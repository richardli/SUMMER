# Construct the scaled Type IV precision and identifying constraints used by
# smoothDirect() and smoothSurvey(). The latent vector is ordered by region,
# with all time points for a region stored consecutively.
.type4_inla_components <- function(Amat, n.time, rw.order = 1L,
                                   constr = TRUE) {
  W <- as.matrix(Amat)
  if (nrow(W) != ncol(W) || anyNA(W) || any(!is.finite(W))) {
    stop("Amat must be a finite square adjacency matrix.", call. = FALSE)
  }
  diag(W) <- 0
  W <- (W != 0) * 1
  if (!isTRUE(all.equal(W, t(W), tolerance = 0))) {
    stop(
      "Amat must have symmetric neighbour links for a Besag interaction.",
      call. = FALSE
    )
  }

  n.region <- nrow(W)
  Q.space <- diag(rowSums(W)) - W
  inla.rw <- utils::getFromNamespace("inla.rw", "INLA")
  inla.scale.model.bym <- utils::getFromNamespace(
    "inla.scale.model.bym", "INLA"
  )
  inla.bym.constr.internal <- utils::getFromNamespace(
    "inla.bym.constr.internal", "INLA"
  )

  component.info <- inla.bym.constr.internal(
    Q.space, adjust.for.con.comp = TRUE
  )
  Q.space <- inla.scale.model.bym(Q.space, adjust.for.con.comp = TRUE)
  Q.time <- inla.rw(
    n.time, order = rw.order, scale.model = TRUE, sparse = TRUE
  )

  extra <- NULL
  reference <- integer()
  if (isTRUE(constr)) {
    component.A <- as.matrix(component.info$constr$A)
    if (nrow(component.A)) {
      reference <- apply(component.A, 1L, function(x) which(x != 0)[1L])
    }

    temporal.nodes <- setdiff(seq_len(n.region), reference)
    temporal.A <- matrix(0, length(temporal.nodes), n.time * n.region)
    for (i in seq_along(temporal.nodes)) {
      r <- temporal.nodes[i]
      temporal.A[i, ((r - 1L) * n.time + 1L):(r * n.time)] <- 1
    }

    spatial.A <- matrix(
      0, n.time * nrow(component.A), n.time * n.region
    )
    if (nrow(component.A)) {
      for (cc in seq_len(nrow(component.A))) {
        nodes <- which(component.A[cc, ] != 0)
        for (tt in seq_len(n.time)) {
          cols <- (nodes - 1L) * n.time + tt
          spatial.A[tt + (cc - 1L) * n.time, cols] <- 1
        }
      }
    }

    A <- rbind(temporal.A, spatial.A)
    if (nrow(A)) extra <- list(A = A, e = rep(0, nrow(A)))
  }

  list(
    Q = Q.space %x% Q.time,
    extraconstr = extra,
    rankdef = n.time * n.region -
      (n.time - rw.order) * (n.region - component.info$rankdef),
    component.info = component.info,
    reference = reference
  )
}
