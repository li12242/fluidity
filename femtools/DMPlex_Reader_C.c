#include "DMPlex_Reader_C.h"

PetscErrorCode dmplex_get_gmsh_plex(const char filename[], DM *plex)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMPlexCreateGmshFromFile(PETSC_COMM_WORLD, filename, PETSC_TRUE, plex);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode dmplex_get_mesh_connectivity(DM plex, PetscInt nnodes, PetscInt loc, PetscInt *ndglno)
{
  PetscInt c, cStart, cEnd, vStart, vEnd, idx, ci, nclosure, point;
  PetscInt *closure=NULL;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMPlexGetHeightStratum(plex, 0, &cStart, &cEnd);CHKERRQ(ierr);
  ierr = DMPlexGetDepthStratum(plex, 0, &vStart, &vEnd);CHKERRQ(ierr);
  for (idx=0, c=cStart; c<cEnd; c++) {
    ierr = DMPlexGetTransitiveClosure(plex, c, PETSC_TRUE, &nclosure, &closure);CHKERRQ(ierr);
    for (ci=0; ci<nclosure; ci++) {
      point = closure[2*ci];
      if (vStart <= point && point < vEnd) ndglno[idx++] = point - vStart + 1;
    }
  }
  if (closure) {
    ierr = DMPlexRestoreTransitiveClosure(plex, c, PETSC_TRUE, &nclosure, &closure);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode dmplex_get_num_surface_facets(DM plex, const char label_name[], PetscInt *nfacets)
{
  PetscInt v, nvalues, npoints;
  const PetscInt *values;
  DMLabel label;
  IS valueIS;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMPlexGetLabel(plex, label_name, &label);CHKERRQ(ierr);
  if (!label) {*nfacets = 0; PetscFunctionReturn(0);}
  ierr = DMLabelGetNumValues(label, &nvalues);CHKERRQ(ierr);
  ierr = DMLabelGetValueIS(label, &valueIS);CHKERRQ(ierr);
  ierr = ISGetIndices(valueIS, &values);CHKERRQ(ierr);
  for (*nfacets=0, v=0; v<nvalues; v++) {
    ierr = DMLabelGetStratumSize(label, values[v], &npoints);CHKERRQ(ierr);
    *nfacets += npoints;
  }
  ierr = ISRestoreIndices(valueIS, &values);CHKERRQ(ierr);
  ierr = ISDestroy(&valueIS);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode dmplex_get_surface_connectivity(DM plex, const char label_name[], PetscInt nfacets, PetscInt sloc, PetscInt *sndglno, PetscInt *boundary_ids)
{
  PetscInt        v, vStart, vEnd, nvalues, p, npoints, idx, ci, nclosure, vertex, nvertices;
  const PetscInt *values, *points;
  DMLabel         label;
  IS              valueIS, pointIS;
  PetscInt       *closure = NULL;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  ierr = DMPlexGetDepthStratum(plex, 0, &vStart, &vEnd);CHKERRQ(ierr);
  ierr = DMPlexGetLabel(plex, label_name, &label);CHKERRQ(ierr);
  if (!label) PetscFunctionReturn(0);
  ierr = DMLabelGetNumValues(label, &nvalues);CHKERRQ(ierr);
  ierr = DMLabelGetValueIS(label, &valueIS);CHKERRQ(ierr);
  ierr = ISGetIndices(valueIS, &values);CHKERRQ(ierr);
  /* Loop over all marker values in the supplied label */
  for (idx = 0 , v = 0; v < nvalues; v++) {
    ierr = DMLabelGetStratumSize(label, values[v], &npoints);CHKERRQ(ierr);
    ierr = DMLabelGetStratumIS(label, values[v], &pointIS);CHKERRQ(ierr);
    ierr = ISGetIndices(pointIS, &points);CHKERRQ(ierr);
    for (p = 0; p < npoints; p++) {
      /* Derive vertices for each marked facet */
      ierr = DMPlexGetTransitiveClosure(plex, points[p], PETSC_TRUE, &nclosure, &closure);CHKERRQ(ierr);
      for (nvertices = 0, ci = 0; ci < nclosure; ci++) {
        vertex = closure[2*ci];
        if (vStart <= vertex && vertex < vEnd) {
          sndglno[idx*sloc+nvertices++] = vertex - vStart + 1;
        }
      }
      /* Store associated boundary ID */
      boundary_ids[idx] = values[v];
      idx++;
    }
    ierr = ISRestoreIndices(pointIS, &points);CHKERRQ(ierr);
    ierr = ISDestroy(&pointIS);CHKERRQ(ierr);
  }
  ierr = ISRestoreIndices(valueIS, &values);CHKERRQ(ierr);
  ierr = ISDestroy(&valueIS);CHKERRQ(ierr);
  if (closure) {
    ierr = DMPlexRestoreTransitiveClosure(plex, 0, PETSC_TRUE, &nclosure, &closure);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode dmplex_get_halo_send_recv(DM plex, PetscSF *sf, PetscInt height,
                                         PetscInt nprocs, PetscInt *nowned,
                                         PetscInt *nsend, IS *sendIS,
                                         PetscInt *nrecv, IS *recvIS) {
  MPI_Comm           comm;
  PetscInt           p, vStart, vEnd, r, nranks, rstart, nleaves, recvSize, sendSize, offset;
  PetscInt          *pidx, *recv_arr, *send_arr;
  const PetscInt    *ilocal, *rranks;
  const PetscSFNode *iremote;
  PetscSection       recvSection, sendSection, rootSection, leafSection;
  IS                 rootRanks, leafRanks;
  PetscErrorCode     ierr;

  PetscFunctionBegin;
  ierr = PetscObjectGetComm((PetscObject) plex, &comm);CHKERRQ(ierr);
  ierr = DMPlexGetDepthStratum(plex, height, &vStart, &vEnd);CHKERRQ(ierr);
  *nowned = vEnd - vStart;

  /* Derive receiver mapping from point SF */
  ierr = PetscCalloc1(nprocs, &pidx);CHKERRQ(ierr);
  for (p = 0; p < nprocs; p++) {nrecv[p] = 0;}
  ierr = PetscSFGetGraph(*sf, NULL, &nleaves, &ilocal, &iremote);CHKERRQ(ierr);
  ierr = PetscSectionCreate(comm, &recvSection);CHKERRQ(ierr);
  ierr = PetscSectionSetChart(recvSection, 0, nprocs);CHKERRQ(ierr);
  for (p = 0; p < nleaves; p++) {
    if (vStart <= ilocal[p] && ilocal[p] < vEnd) {
      ierr = PetscSectionAddDof(recvSection, iremote[p].rank, 1);CHKERRQ(ierr);
      (*nowned)--;
    }
  }
  ierr = PetscSectionSetUp(recvSection);CHKERRQ(ierr);
  /* */
  ierr = PetscSectionGetStorageSize(recvSection, &recvSize);CHKERRQ(ierr);
  ierr = PetscMalloc1(recvSize, &recv_arr);CHKERRQ(ierr);
  for (p = 0; p < nleaves; p++) {
    if (vStart <= ilocal[p] && ilocal[p] < vEnd) {
      const PetscInt proc = iremote[p].rank;
      ierr = PetscSectionGetOffset(recvSection, proc, &offset);CHKERRQ(ierr);
      recv_arr[offset+pidx[proc]++] = ilocal[p] - vStart + 1;
    }
  }
  for (p = 0; p < nprocs; p++) {
    ierr = PetscSectionGetDof(recvSection, p, &(nrecv[p]));CHKERRQ(ierr);
  }
  ierr = ISCreateGeneral(comm, recvSize, recv_arr, PETSC_OWN_POINTER, recvIS);CHKERRQ(ierr);

  /* Derive sender mapping from root/rank mapping */
  PetscSF sfPoint;
  ierr = DMGetPointSF(plex, &sfPoint);CHKERRQ(ierr);
  ierr = PetscObjectReference((PetscObject)sfPoint);CHKERRQ(ierr);
  ierr = DMSetPointSF(plex, *sf);CHKERRQ(ierr);

  ierr = PetscSectionCreate(comm, &rootSection);CHKERRQ(ierr);
  ierr = PetscSectionCreate(comm, &leafSection);CHKERRQ(ierr);
  ierr = DMPlexDistributeOwnership(plex, rootSection, &rootRanks, leafSection, &leafRanks);CHKERRQ(ierr);
  ierr = ISGetIndices(rootRanks, &rranks);CHKERRQ(ierr);
  ierr = DMSetPointSF(plex, sfPoint);CHKERRQ(ierr);
  ierr = PetscSFDestroy(&sfPoint);CHKERRQ(ierr);

  /* */
  ierr = PetscSectionCreate(comm, &sendSection);CHKERRQ(ierr);
  ierr = PetscSectionSetChart(sendSection, 0, nprocs);CHKERRQ(ierr);
  for (p = vStart; p < vEnd; p++) {
    ierr = PetscSectionGetDof(rootSection, p, &nranks);CHKERRQ(ierr);
    ierr = PetscSectionGetOffset(rootSection, p, &rstart);CHKERRQ(ierr);
    for (r = 0; r < nranks; r++) {
      ierr = PetscSectionAddDof(sendSection, rranks[rstart+r], 1);CHKERRQ(ierr);
    }
  }
  ierr = PetscSectionSetUp(sendSection);CHKERRQ(ierr);
  ierr = PetscSectionGetStorageSize(sendSection, &sendSize);CHKERRQ(ierr);
  /* */
  for (p = 0; p < nprocs; p++) {nsend[p] = 0; pidx[p] = 0;}
  ierr = PetscMalloc1(sendSize, &send_arr);CHKERRQ(ierr);
  for (p = vStart; p < vEnd; p++) {
    ierr = PetscSectionGetDof(rootSection, p, &nranks);CHKERRQ(ierr);
    ierr = PetscSectionGetOffset(rootSection, p, &rstart);CHKERRQ(ierr);
    for (r = 0; r < nranks; r++) {
      const PetscInt proc = rranks[rstart+r];
      ierr = PetscSectionGetOffset(sendSection, proc, &offset);CHKERRQ(ierr);
      send_arr[offset+pidx[proc]++] = p - vStart + 1;
    }
  }
  for (p = 0; p < nprocs; p++) {
    ierr = PetscSectionGetDof(sendSection, p, &(nsend[p]));CHKERRQ(ierr);
  }
  ierr = ISCreateGeneral(comm, sendSize, send_arr, PETSC_OWN_POINTER, sendIS);CHKERRQ(ierr);

  ierr = PetscSectionDestroy(&recvSection);CHKERRQ(ierr);
  ierr = PetscSectionDestroy(&sendSection);CHKERRQ(ierr);
  ierr = PetscSectionDestroy(&rootSection);CHKERRQ(ierr);
  ierr = ISRestoreIndices(rootRanks, &rranks);CHKERRQ(ierr);
  ierr = ISDestroy(&rootRanks);CHKERRQ(ierr);
  ierr = PetscSectionDestroy(&leafSection);CHKERRQ(ierr);
  ierr = ISDestroy(&leafRanks);CHKERRQ(ierr);
  ierr = PetscFree(pidx);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
