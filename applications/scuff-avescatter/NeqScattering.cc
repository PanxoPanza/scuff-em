/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of SCUFF-EM.
 *
 * SCUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * SCUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * NeqScattering.cc  -- Compute non-equilibrium scattering calculations
 *                   -- using Energy/Momentum Transfer PFT method (EMTPFT) 
 *                   -- for a single frequency
 *
 * Francisco Ramirez  -- 01/2020
 *
 */

#include "scuff-avescatter.h"
#include "libscuffInternals.h"
#ifdef USE_OPENMP
#  include <omp.h>
#endif

#define II cdouble(0.0,1.0)
// constant to define exterior index and neq-scattering calculations (F. Ramirez 2020)

#define NUMPFTQ 17    // 2P + 3F + 6T1 + 6T2
#define NUMF    3
#define NUMPFTIS 3*NUMPFTQ
#define NO_SCUFF_TRANSFORM false

/***************************************************************/
/***************************************************************/
/***************************************************************/
void UndoSCUFFMatrixTransformation(HMatrix *M)
{ 
  for (int nr=0; nr<M->NR; nr+=2)
   for (int nc=0; nc<M->NC; nc+=2)
    { M->SetEntry(nr,   nc,   ZVAC*M->GetEntry(nr,   nc)   );
      M->SetEntry(nr,   nc+1, -1.0*M->GetEntry(nr,   nc+1) );
      M->SetEntry(nr+1, nc+1, -1.0*M->GetEntry(nr+1, nc+1)/ZVAC );
    };
}

/***************************************************************/
/* compute the trace formula for the product between two matrix*/
/*		- Gdest is the destiny matrix                          */
/*		- DRMatrix is the Dressed Rytov Matrix (4/PI Factor)   */
/*                                                             */
/*  Note: Gdest matrix comes with scuff transformation         */
/* Francisco Ramirez 2019									   */
/***************************************************************/
double MatrixTraceFlux(HMatrix *GDest, HMatrix* DRMatrix, bool scuff_transform = true)
{
	double scuff_scale[2][2];

	if (scuff_transform) {
		scuff_scale[0][0] = ZVAC;
		scuff_scale[0][1] = -1.0;
		scuff_scale[1][0] = +1.0;
		scuff_scale[1][1] = -1.0 / ZVAC;
	}
	else {
		scuff_scale[0][0] = 1.0;
		scuff_scale[0][1] = 1.0;
		scuff_scale[1][0] = 1.0;
		scuff_scale[1][1] = 1.0;
	}

	// Extract characteristic parameters from the destiny surface
	int DimD = GDest->NR;
	int DimS = DRMatrix->NC;

	double FMPTrace = 0.0; //'two-matrix-product trace'
		// Get trace
	for (int nr = 0; nr<DimD; nr += 2)
		for (int nc = 0; nc<DimS; nc += 2)
			for (int i = 0; i < 2; i++)
				for (int j = 0; j < 2; j++)
					FMPTrace += real(scuff_scale[i][j]
						* GDest->GetEntry(nr + i, nc + j)
						*DRMatrix->GetEntry(nc + j, nr + i));
	return FMPTrace;
}

/****************************************************************/
/*               Compute scattering trace flux					*/
/****************************************************************/
double ScaTraceFlux(SNEQData *SNEQD, HMatrix *GDest, HMatrix* DRMatrix,
	int DestSurface, bool scuff_transform = true)
{
	RWGGeometry *G = SNEQD->G;
	int DimS;
	int OffsetS;
	int NS = G->NumSurfaces;
	int RegionIdx = 0;
	HMatrix *DRblock = 0, *GDblock = 0;
	
	//Log("Computing scattering at region %i",RegionIdx);
	
	if (G->Surfaces[DestSurface]->RegionIndices[0] != RegionIdx &&
		G->Surfaces[DestSurface]->RegionIndices[1] != RegionIdx)
		return 0.0;

	// Extract characteristic parameters from the destiny surface
	int DimD = G->Surfaces[DestSurface]->NumBFs;
	int OffsetD = G->BFIndexOffset[DestSurface];

	double FMPTrace = 0.0; //'two-matrix-product trace'
	for (int nss = 0; nss < NS; nss++) {
		
		// Prepare matrix block
		DimS = G->Surfaces[nss]->NumBFs;
		OffsetS = G->BFIndexOffset[nss];

		if (G->Surfaces[nss]->RegionIndices[0] != RegionIdx &&
		    G->Surfaces[nss]->RegionIndices[1] != RegionIdx)
			continue;

		DRblock = new HMatrix(DimS, DimD, LHM_COMPLEX);
		DRMatrix->ExtractBlock(OffsetS, OffsetD, DRblock);

		GDblock = new HMatrix(DimD, DimS, LHM_COMPLEX);
		GDest->ExtractBlock(OffsetD, OffsetS, GDblock);

		// Get trace
		//Log("...Trace, Region=%i, obj = {%i,%i}",RegionIdx,DestSurface+1,nss+1);
		FMPTrace += MatrixTraceFlux(GDblock, DRblock, scuff_transform);
		
		delete GDblock; GDblock = 0;
		delete DRblock; DRblock = 0;
	}
	return -(1.0 / 2.0) * FMPTrace;
}

/****************************************************************/
/*               Compute absorption trace flux					*/
/****************************************************************/
double AbsTraceFlux(SNEQData *SNEQD, HMatrix* DRMatrix,
	int RegionIdx, bool scuff_transform = true)
{
	RWGGeometry *G = SNEQD->G;
	HMatrix **TExt = SNEQD->TExt;
	HMatrix **TInt = SNEQD->TInt;
	HMatrix **U = SNEQD->U;
	int DimS, DimD;
	int OffsetS, OffsetD;
	int NS = G->NumSurfaces;
	HMatrix *GDblock = 0;
	HMatrix *DRblock = 0;

	double FMPTrace = 0.0; // two-matrix-product trace
	//Log("Computing absorption at Region %i",RegionIdx);
	for (int nb = 0, nsr = 0; nsr < NS; nsr++) {
		// Prepare matrix block
		DimD = G->Surfaces[nsr]->NumBFs;
		OffsetD = G->BFIndexOffset[nsr];

		//---------------------------------------------------------------------
		// Compute elements in the diagonal
		//---------------------------------------------------------------------
		DimS = G->Surfaces[nsr]->NumBFs;
		OffsetS = G->BFIndexOffset[nsr];

		if      (G->Surfaces[nsr]->RegionIndices[0] == RegionIdx)
			GDblock = CopyHMatrixUnpacked(TExt[nsr]);
		else if (G->Surfaces[nsr]->RegionIndices[1] == RegionIdx)
			GDblock = CopyHMatrixUnpacked(TInt[nsr]);

		if (GDblock != 0) {
			DRblock = new HMatrix(DimS, DimD, LHM_COMPLEX);
			DRMatrix->ExtractBlock(OffsetS, OffsetD, DRblock);
			
			//Log("...Trace, Region %i, obj = {%i,%i}",RegionIdx,nsr+1,nsr+1);
			FMPTrace += MatrixTraceFlux(GDblock, DRblock, scuff_transform);
			delete DRblock;
			delete GDblock; GDblock = 0;
		}
		else {// surface nsr is not connected with RegionIdx, skip row
			nb = 0 ? -1 : nb;
			nb += NS - (nsr + 1);
			continue;
		}

		//---------------------------------------------------------------------
		// Compute elements off the diagonal
		//---------------------------------------------------------------------
		for (int nsc = nsr + 1; nsc < NS; nsc++, nb++) {
			//Log("nb = %i", nb);
			DimS = G->Surfaces[nsc]->NumBFs;
			OffsetS = G->BFIndexOffset[nsc];

			//Log("obj%i, ExtRegion %i", nsc + 1, G->Surfaces[nsc]->RegionIndices[0]);
			//Log("obj%i, IntRegion %i", nsc + 1, G->Surfaces[nsc]->RegionIndices[1]);
			if (G->Surfaces[nsc]->RegionIndices[0] == RegionIdx ||
				G->Surfaces[nsc]->RegionIndices[1] == RegionIdx)
				GDblock = CopyHMatrixUnpacked(U[nb]);

			if (GDblock != 0) {
				// flux at upper triangular part
				DRblock = new HMatrix(DimS, DimD, LHM_COMPLEX);
				DRMatrix->ExtractBlock(OffsetS, OffsetD, DRblock);

				//Log("...Trace, Region %i, obj = {%i,%i}",RegionIdx,nsr+1,nsc+1);
				//Log("GDblock, {NR, NC} = {%i, %i}", GDblock->NR, GDblock->NC);
				//Log("DRblock, {NR, NC} = {%i, %i}", DRblock->NR, DRblock->NC);
				FMPTrace += MatrixTraceFlux(GDblock, DRblock, scuff_transform);
				delete DRblock;

				// flux at lower triangular part
				GDblock->Transpose();
				DRblock = new HMatrix(DimD, DimS, LHM_COMPLEX);
				DRMatrix->ExtractBlock(OffsetD, OffsetS, DRblock);
				
				//Log("...Trace, Region %i, obj = {%i,%i}",RegionIdx,nsc+1,nsr+1);
				FMPTrace += MatrixTraceFlux(GDblock, DRblock, scuff_transform);
				delete DRblock;
				delete GDblock;	GDblock = 0;
			}
		}
	}
	return -(1.0 / 2.0) * FMPTrace;
}

double GetTraceFlux(SNEQData *SNEQD, HMatrix *GDest, HMatrix* DRMatrix,
	int DestSurface, int RegionIdx = 0, bool scuff_transform = true)
{
	if (RegionIdx == 0)
		return ScaTraceFlux(SNEQD, GDest, DRMatrix, DestSurface, scuff_transform);
	else
		return AbsTraceFlux(SNEQD, DRMatrix, RegionIdx, scuff_transform);

	return 0.0;
}


/***************************************************************/
/* Compute Anti-Hermitian part of a matrix dT0, whose          */ 
/* elements are anti-symmetrical.                              */
/*                   - i*asymA = (A - A*)/2i                   */
/* F. Ramirez 2019											   */
/***************************************************************/
void iAntiHermitianPart(RWGGeometry *G, HMatrix *A, HMatrix *asymA, 
	cdouble SCALE_FACT = 1.0)
{
	int TotalEdges = G->TotalEdges;
	cdouble asymQKK_ab, asymQKN_ab, asymQNK_ab, asymQNN_ab;
	cdouble asymQKK_ba, asymQKN_ba, asymQNK_ba, asymQNN_ba;

	asymA->Zero();
	for (int neaTot = 0; neaTot < TotalEdges; neaTot++)
		for (int nebTot = neaTot; nebTot < TotalEdges; nebTot++) {

			// Retrieve index from surface A
			int nsa, nea, KNIndexA;
			G->ResolveEdge(neaTot, &nsa, &nea, &KNIndexA);

			// Retrieve index from surface B
			int nsb, neb, KNIndexB;
			G->ResolveEdge(nebTot, &nsb, &neb, &KNIndexB);

			// Extract QF elements from upper triangular part
			cdouble QKK_ab = ZVAC * A->GetEntry(2 * neaTot + 0, 2 * nebTot + 0);
			cdouble QKN_ab = - 1.0 * A->GetEntry(2 * neaTot + 0, 2 * nebTot + 1);
			cdouble QNK_ab = A->GetEntry(2 * neaTot + 1, 2 * nebTot + 0);
			cdouble QNN_ab = - 1.0 * A->GetEntry(2 * neaTot + 1, 2 * nebTot + 1) / ZVAC;

			// Extract QF elements from lower triangular part
			cdouble QKK_ba = ZVAC * A->GetEntry(2 * nebTot + 0, 2 * neaTot + 0);
			cdouble QKN_ba = - 1.0 * A->GetEntry(2 * nebTot + 0, 2 * neaTot + 1);
			cdouble QNK_ba = A->GetEntry(2 * nebTot + 1, 2 * neaTot + 0);
			cdouble QNN_ba = - 1.0 * A->GetEntry(2 * nebTot + 1, 2 * neaTot + 1) / ZVAC;

			// Scale elements by ALPHA
			QKK_ab = SCALE_FACT*QKK_ab; QKK_ba = SCALE_FACT*QKK_ba;
			QKN_ab = SCALE_FACT*QKN_ab; QKN_ba = SCALE_FACT*QKN_ba;
			QNK_ab = SCALE_FACT*QNK_ab; QNK_ba = SCALE_FACT*QNK_ba;
			QNN_ab = SCALE_FACT*QNN_ab; QNN_ba = SCALE_FACT*QNN_ba;

			if (nsa == nsb) { // same surface
				if (nea == neb) { // same edges
								  // Upper/Lower triangular part, ab elements
					asymQKK_ab = 0.0 + II*0.0; asymQKK_ba = 0.0 + II*0.0;
					asymQKN_ab = 0.0 + II*0.0; asymQKN_ba = 0.0 + II*0.0;
					asymQNK_ab = 0.0 + II*0.0; asymQNK_ba = 0.0 + II*0.0;
					asymQNN_ab = 0.0 + II*0.0; asymQNN_ba = 0.0 + II*0.0;

				}
				else {  // different edges, same surface
						// Upper triangular part, ab elements
					asymQKK_ab = -II*real(QKK_ab);
					asymQKN_ab = +imag(QKN_ab);
					asymQNK_ab = +imag(QNK_ab);
					asymQNN_ab = -II*real(QNN_ab);

					// Lower triangular part, ba elements
					asymQKK_ba = -II*real(QKK_ba);
					asymQKN_ba = +imag(QKN_ba);
					asymQNK_ba = +imag(QNK_ba);
					asymQNN_ba = -II*real(QNN_ba);
				}
			}
			else { // different surfaces
				   // Upper triangular part, ab elements
				asymQKK_ab = -0.5*II*(QKK_ab - conj(QKK_ba));
				asymQKN_ab = -0.5*II*(QKN_ab - conj(QNK_ba));
				asymQNK_ab = -0.5*II*(QNK_ab - conj(QKN_ba));
				asymQNN_ab = -0.5*II*(QNN_ab - conj(QNN_ba));

				// Lower triangular part, ba elements
				asymQKK_ba = -0.5*II*(QKK_ba - conj(QKK_ab));
				asymQKN_ba = -0.5*II*(QKN_ba - conj(QNK_ab));
				asymQNK_ba = -0.5*II*(QNK_ba - conj(QKN_ab));
				asymQNN_ba = -0.5*II*(QNN_ba - conj(QNN_ab));
			}

			/*..................................................*/
			/* Assemble asym(Qforce)                            */
			/*..................................................*/
			// Upper triangular part, ab elements
			asymA->SetEntry(2 * neaTot + 0, 2 * nebTot + 0, asymQKK_ab);
			asymA->SetEntry(2 * neaTot + 0, 2 * nebTot + 1, asymQKN_ab);
			asymA->SetEntry(2 * neaTot + 1, 2 * nebTot + 0, asymQNK_ab);
			asymA->SetEntry(2 * neaTot + 1, 2 * nebTot + 1, asymQNN_ab);

			// Lower triangular part, ba elements
			asymA->SetEntry(2 * nebTot + 0, 2 * neaTot + 0, asymQKK_ba);
			asymA->SetEntry(2 * nebTot + 0, 2 * neaTot + 1, asymQKN_ba);
			asymA->SetEntry(2 * nebTot + 1, 2 * neaTot + 0, asymQNK_ba);
			asymA->SetEntry(2 * nebTot + 1, 2 * neaTot + 1, asymQNN_ba);
		}
}

/***************************************************************/
/* Compute Hermitian part of a matrix A						   */
/*                   symA = (A + A*)/2                         */
/* F. Ramirez 2019											   */
/***************************************************************/
void HermitianPart(HMatrix *A, HMatrix *symA, int OffsetS, 
	cdouble SCALE_FACT = 1.0) {

	int NBFSr = A->NR;
	int NBFSc = A->NC;
	/***************************************************************/
	/* stamp Sym(T_s) = (T_s + T_s^\dagger) / 2                    */
	/* into the sth diagonal block of DRMatrix,                    */
	/* undoing the SCUFF matrix transformation along the way.      */
	/***************************************************************/
	symA->Zero();
	for (int a = 0; a < (NBFSr / 2); a++)
		for (int b = a; b < (NBFSc / 2); b++)
		{
			cdouble TEEab = ZVAC * A->GetEntry(2 * a + 0, 2 * b + 0);
			cdouble TEMab = -1.0 * A->GetEntry(2 * a + 0, 2 * b + 1);
			cdouble TMEab = A->GetEntry(2 * a + 1, 2 * b + 0);
			cdouble TMMab = -1.0 * A->GetEntry(2 * a + 1, 2 * b + 1) / ZVAC;

			cdouble TEEba = ZVAC * A->GetEntry(2 * b + 0, 2 * a + 0);
			cdouble TEMba = -1.0 * A->GetEntry(2 * b + 0, 2 * a + 1);
			cdouble TMEba = A->GetEntry(2 * b + 1, 2 * a + 0);
			cdouble TMMba = -1.0 * A->GetEntry(2 * b + 1, 2 * a + 1) / ZVAC;

			cdouble SymTEE = 0.5*SCALE_FACT*(TEEab + conj(TEEba));
			cdouble SymTEM = 0.5*SCALE_FACT*(TEMab + conj(TMEba));
			cdouble SymTME = 0.5*SCALE_FACT*(TMEab + conj(TEMba));
			cdouble SymTMM = 0.5*SCALE_FACT*(TMMab + conj(TMMba));

			symA->SetEntry(OffsetS + 2 * a + 0, OffsetS + 2 * b + 0, SymTEE);
			symA->SetEntry(OffsetS + 2 * b + 0, OffsetS + 2 * a + 0, conj(SymTEE));

			symA->SetEntry(OffsetS + 2 * a + 0, OffsetS + 2 * b + 1, SymTEM);
			symA->SetEntry(OffsetS + 2 * b + 1, OffsetS + 2 * a + 0, conj(SymTEM));

			symA->SetEntry(OffsetS + 2 * a + 1, OffsetS + 2 * b + 0, SymTME);
			symA->SetEntry(OffsetS + 2 * b + 0, OffsetS + 2 * a + 1, conj(SymTME));

			symA->SetEntry(OffsetS + 2 * a + 1, OffsetS + 2 * b + 1, SymTMM);
			symA->SetEntry(OffsetS + 2 * b + 1, OffsetS + 2 * a + 1, conj(SymTMM));
		};
}

/***************************************************************/
/* Compute the dressed Rytov matrix for sources contained in   */
/* SourceSurface. The matrix is stored in the DRMatrix         */
/* field of the SNEQD structure.                               */
/***************************************************************/
#define RYTOVPF (-4.0/M_PI)
void AdvancedDRMatrix(SNEQData *SNEQD, int SourceSurface)
{
	RWGGeometry *G = SNEQD->G;
	HMatrix *M = SNEQD->M;
	HMatrix *DR = SNEQD->DRMatrix;
	int OffsetS = 0;

	if (SourceSurface >= 0) { // compute DR matrix for interior surfaces (F. Ramirez 2019)
		Log("computing DR matrix for internal object %i", SourceSurface + 1);
		OffsetS = G->BFIndexOffset[SourceSurface];
		HermitianPart(SNEQD->TInt[SourceSurface], DR, OffsetS, RYTOVPF);
	}
	else if (SourceSurface == AV_SCATTERING) { // DR matrix for exterior (F. Ramirez 2019)
		Log("computing DR matrix for exterior (F.Ramirez 2019)...");
		HermitianPart(SNEQD->T0, DR, OffsetS, RYTOVPF);
	}
	else if (SourceSurface == AV_ASYM_X) { //DR matrix for gcos_x (F. Ramirez 2019)
		Log("computing DR matrix for gcos-x (F.Ramirez 2019)...");
		iAntiHermitianPart(G, SNEQD->T0, DR, RYTOVPF);
	}
	else if (SourceSurface == AV_ASYM_Y) { // DR matrix for gcos_y (F. Ramirez 2019)
		Log("computing DR matrix for gcos-y (F.Ramirez 2019)...");
		iAntiHermitianPart(G, SNEQD->T0, DR, RYTOVPF);
	}
	else if (SourceSurface == AV_ASYM_Z) { // DR matrix for gcos_z (F. Ramirez 2019)
		Log("computing DR matrix for gcos-z (F.Ramirez 2019)...");
		iAntiHermitianPart(G, SNEQD->T0, DR, RYTOVPF);
	}

	/***************************************************************/
	/* set DR = W * DR * W' ****************************************/
	/* by computing DR = M \ (M \ DR)'       ***********************/
	/***************************************************************/
	M->LUSolve(DR); // M \ DR        = W * DR
	DR->Adjoint();  // (M \ DR)'     = DR * W'
	M->LUSolve(DR); // M \ (M \ DR)' = W * DR * W' 
	Log("...done with DR matrix");
}

// Compute Green function gradients from EMT method (F. Ramirez 2020)
void GetForceIntegrals(RWGGeometry *G,
	int nsa, int nea, int nsb, int neb,
	cdouble Omega, cdouble k,
	cdouble EpsR, cdouble MuR,
	cdouble PFTIs[3 * NUMPFTQ])
{
	RWGSurface *SA = G->Surfaces[nsa];
	RWGSurface *SB = G->Surfaces[nsb];
	double rRel;
	AssessBFPair(SA, nea, SB, neb, &rRel);

	GetGCMEArgStruct MyArgs, *Args = &MyArgs;
	InitGetGCMEArgs(Args);
	Args->nsa = nsa;
	Args->nsb = nsb;
	Args->NumRegions = 1;
	Args->k[0] = k;
	Args->NeedForce = true;
	Args->FIBBICache = (nsa == nsb) ? G->FIBBICaches[nsa] : 0;
	cdouble GabArray[2][NUMGCMES];
	cdouble ikCabArray[2][NUMGCMES];
	GetGCMatrixElements(G, Args, nea, neb, GabArray, ikCabArray);

	cdouble *QKK = PFTIs + 0 * NUMPFTQ;
	cdouble *QNN = PFTIs + 1 * NUMPFTQ;
	cdouble *QKNmNK = PFTIs + 2 * NUMPFTQ;
	double FTPreFac = TENTHIRDS;
	for (int Mu = 0; Mu<NUMF; Mu++)
	{
		QKK[PFT_XFORCE + Mu] = FTPreFac*MuR*GabArray[0][GCME_FX + Mu];
		QNN[PFT_XFORCE + Mu] = FTPreFac*EpsR*GabArray[0][GCME_FX + Mu];
		QKNmNK[PFT_XFORCE + Mu] = FTPreFac*ikCabArray[0][GCME_FX + Mu] / Omega;
	}
}

/***************************************************************/
/* Store elements of KNmatrix and QForce		               */
/* F. Ramirez 2019									           */
/***************************************************************/
void StoreMatrix(RWGGeometry *G, int nsa, int nea, int nsb, int neb,
	cdouble Omega, int RegionIndex, bool UseSymmetry, HMatrix *QForce,
	double Sign, int FORCE_IDX)
{
	cdouble EpsR, MuR;
	G->RegionMPs[RegionIndex]->GetEpsMu(Omega, &EpsR, &MuR);
	cdouble k = Omega*sqrt(EpsR*MuR);
	cdouble PFTIs[3* NUMPFTQ];
	cdouble *QKK = PFTIs + 0 * NUMPFTQ;
	cdouble *QNN = PFTIs + 1 * NUMPFTQ;
	cdouble *QKNmNK = PFTIs + 2 * NUMPFTQ;
	GetForceIntegrals(G, nsa, nea, nsb, neb,
		Omega, k, EpsR, MuR, PFTIs);

	
	int OffsetSa = G->BFIndexOffset[nsa];
	int OffsetSb = G->BFIndexOffset[nsb];
	
	/*..........................................................*/
	/*........ Save Qforce elements (F. Ramirez 2019) ..........*/
	/*..........................................................*/
	int a_index = OffsetSa + 2 * nea;
	int b_intex = OffsetSb + 2 * neb;

	//for (int nF = 0; nF < NUMF; nF++) {
		QForce->SetEntry(a_index + 0, b_intex + 0, +Sign*II*QKK[FORCE_IDX]);
		QForce->SetEntry(a_index + 0, b_intex + 1, -Sign*QKNmNK[FORCE_IDX]);
		QForce->SetEntry(a_index + 1, b_intex + 0, -Sign*QKNmNK[FORCE_IDX]);
		QForce->SetEntry(a_index + 1, b_intex + 1, -Sign*II*QNN[FORCE_IDX]);

		/*************************************************************/
		/*   If symmetric problem, fill the lower triangular part    */
		/*************************************************************/
		if (UseSymmetry) {
			QForce->SetEntry(b_intex + 0, a_index + 0, -Sign*II*QKK[FORCE_IDX]);
			QForce->SetEntry(b_intex + 1, a_index + 0, +Sign*QKNmNK[FORCE_IDX]);
			QForce->SetEntry(b_intex + 0, a_index + 1, +Sign*QKNmNK[FORCE_IDX]);
			QForce->SetEntry(b_intex + 1, a_index + 1, +Sign*II*QNN[FORCE_IDX]);
		}
	//}
}

void gradT0_Matrix(RWGGeometry *G, SNEQData *SNEQD, cdouble Omega, HMatrix *dT0, int FORCE_IDX)
{
	PFTOptions *PFTOpts = &(SNEQD->PFTOpts);
	bool Interior = PFTOpts->Interior;
		/*--------------------------------------------------------------*/
		/*- loop over all edge pairs to get scattered PFT contributions */
		/*--------------------------------------------------------------*/
		int NT = 1;
#ifdef USE_OPENMP 
	NT = GetNumThreads();
#endif

	/*--------------------------------------------------------------*/
	/*- multithreaded loop over all basis functions on all surfaces-*/
	/*--------------------------------------------------------------*/
	int TotalEdges = G->TotalEdges;
	bool UseSymmetry = true;
	char *s = getenv("SCUFF_EMTPFT_SYMMETRY");
	if (s && s[0] == '0')
	{
		UseSymmetry = false;
		Log("Not using symmetry in EMTPFT calculation.");
	};
	//TODO insert here a check for any nested objects and automatically 
	//     disable symmetry if present

	int nsaOnly = -1; CheckEnv("SCUFF_EMTPFT_NSAONLY", &nsaOnly);
	int nsbOnly = -1; CheckEnv("SCUFF_EMTPFT_NSBONLY", &nsbOnly);

#ifdef USE_OPENMP
	Log("EMT OpenMP multithreading (%i threads)", NT);
#pragma omp parallel for schedule(dynamic,1), num_threads(NT)
#endif
	for (int neaTot = 0; neaTot<TotalEdges; neaTot++)
		for (int nebTot = (UseSymmetry ? neaTot : 0); nebTot<TotalEdges; nebTot++) {
			if (nebTot == (UseSymmetry ? neaTot : 0))
				LogPercent(neaTot, TotalEdges, 10);

			int nsa, nea, KNIndexA;
			RWGSurface *SA = G->ResolveEdge(neaTot, &nsa, &nea, &KNIndexA);
			int RegionIndex = SA->RegionIndices[Interior ? 1 : 0];
			if (RegionIndex == -1) continue; // no interior PFT for PEC bodies

			int nsb, neb, KNIndexB;
			RWGSurface *SB = G->ResolveEdge(nebTot, &nsb, &neb, &KNIndexB);

			if ((nsaOnly != -1 && nsa != nsaOnly) || (nsbOnly != -1 && nsb != nsbOnly))
				continue;

			if (SA->RegionIndices[0] != 0) // compute only interactions at Exterior
				continue;

			double Sign = 0.0;
			if (nsa == nsb)
				Sign = Interior ? -1.0 : +1.0;
			else if (SA->RegionIndices[0] == SB->RegionIndices[0]) // A, B live in same region
				Sign = (Interior ? 0.0 : 1.0);
			else if (SA->RegionIndices[0] == SB->RegionIndices[1]) // A contained in B
				Sign = Interior ? 0.0 : -1.0;
			else if (SA->RegionIndices[1] == SB->RegionIndices[0]) // B contained in A
				Sign = Interior ? 1.0 : 0.0;

			if (Sign == 0.0) // B does not contribute to PFT on A
				continue;

			/*..........................................................*/
			/*  Store elements of KNmatrix and QForce (F. Ramirez 2019) */
			/*..........................................................*/
			StoreMatrix(G, nsa, nea, nsb, neb, Omega,
				RegionIndex, UseSymmetry, dT0, Sign, FORCE_IDX);
			/*..........................................................*/
		}; // end of multithreaded loop
}

/***************************************************************/
/* Evaluate average absorption and scattering flux             */
/***************************************************************/
void AverageScatteringFlux(SNEQData *SNEQD, cdouble Omega, HMatrix *PFTMatrix)
{
	RWGGeometry *G = SNEQD->G;
	int NS = G->NumSurfaces;
	HMatrix *DRMatrix = SNEQD->DRMatrix;
	HMatrix *T0 = SNEQD->T0;
	//HMatrix *dT0 = SNEQD->dT0;
	HMatrix *asymdG = new HMatrix(G->TotalBFs, G->TotalBFs, LHM_COMPLEX);
	cdouble ScatFlux=0.0, AbsFlux=0.0; 
	cdouble *AsymFlux = new cdouble[NS];
	int RegionIdx;
	//cdouble Aflux; // test
	

	// compute scattering flux for the source (F. Ramirez 2019)
	cdouble EpsR, MuR;
	G->RegionMPs[EXTERIOR_INDEX]->GetEpsMu(Omega, &EpsR, &MuR);
	cdouble khost = Omega*sqrt(EpsR*MuR);
	double kPI2 = real((khost / M_PI)*conj(khost/ M_PI));

	/*----------------------------------------------------------*/
	/*             Compute average sca. & abs. flux             */
	/*----------------------------------------------------------*/
	Log("Computing average scattering and absorption");
	// Compute Dresed Rytov matrix for exterior
	AdvancedDRMatrix(SNEQD, AV_SCATTERING);

	for (int nsd = 0; nsd < NS; nsd++) {
		// Get av. scattering flux
		ScatFlux = GetTraceFlux(SNEQD, T0, DRMatrix, nsd); 
		PFTMatrix->SetEntry(nsd, PFT_PSCAT, ScatFlux / kPI2);

		// Initialize scat<cos> flux
		AsymFlux[nsd] = 0.0;

		// Get av absorption flux
		RegionIdx = G->Surfaces[nsd]->RegionIndices[1];
		if (RegionIdx == 0)
			RegionIdx = G->Surfaces[nsd]->RegionIndices[0];

		AbsFlux = GetTraceFlux(SNEQD, 0, DRMatrix, nsd, RegionIdx);
		PFTMatrix->SetEntry(nsd, PFT_PABS, AbsFlux / kPI2);
		
	}

	/*----------------------------------------------------------*/
	/*             Compute average asymmetry factor             */
	/*----------------------------------------------------------*/
	Log("Computing asymetry factor");
	// Get neq scattering flux
	for (int xyz = 0; xyz < NUM_ASYM; xyz++){
		// compute gradT0
		T0->Zero();
		Log("...Assembling dT0 matrix");
		gradT0_Matrix(G, SNEQD, Omega, T0, PFT_XFORCE + xyz);

		// Compute Dresed Rytov matrix for grad(GC)
		AdvancedDRMatrix(SNEQD, AV_ASYM_X + xyz);
		iAntiHermitianPart(G, T0, asymdG);
		for (int nsd = 0; nsd < NS; nsd++){
			if (G->Surfaces[nsd]->RegionIndices[0] != 0) continue;
			AsymFlux[nsd] += GetTraceFlux(SNEQD, asymdG, DRMatrix, nsd, 0, 
				NO_SCUFF_TRANSFORM);
		}
	}

	// speed of light in medium (x 1E9 m/s)
	cdouble cR = 1.0 / (TENTHIRDS*sqrt(EpsR*MuR));
	for (int nsd = 0; nsd < NS; nsd++) {
	    // Compute asymmetry factor
		AsymFlux[nsd] *= real(cR*conj(cR))/ (kPI2);
		PFTMatrix->SetEntry(nsd, PFT_XFORCE, AsymFlux[nsd]);
	}

	delete[] AsymFlux;
	delete asymdG;
}


/***************************************************************/
/***************************************************************/
/***************************************************************/
void NeqScattering(SNEQData *SNEQD, cdouble Omega, double *kBloch)
{
  /***************************************************************/
  /* extract fields from SNEQData structure **********************/
  /***************************************************************/
  RWGGeometry *G      = SNEQD->G;
  HMatrix *M          = SNEQD->M;
  HMatrix **TExt      = SNEQD->TExt;
  HMatrix **TInt      = SNEQD->TInt; 
  HMatrix **U         = SNEQD->U;
  HMatrix *T0         = SNEQD->T0; // neq sca (F. Ramirez 2019)
  //HMatrix *dT0       = SNEQD->dT0; // neq sca (F. Ramirez 2019)
  int NS              = SNEQD->G->NumSurfaces;
  //char *FileBase      = SNEQD->FileBase;

  Log("Computing neq quantities at omega=%s...",z2s(Omega));

  /***************************************************************/
  /* preinitialize an argument structure for the BEM matrix      */
  /* block assembly routine                                      */
  /***************************************************************/
  GetSSIArgStruct MyGSSIArgStruct, *Args=&MyGSSIArgStruct;
  InitGetSSIArgs(Args);
  Args->G         = G;
  Args->Omega     = Omega;

  /***************************************************************/
  /* before entering the loop over transformations, we first     */
  /* assemble the (transformation-independent) T matrix blocks.  */
  /***************************************************************/
  for(int nr=0; nr<G->NumRegions; nr++)
	  G->RegionMPs[nr]->Zero();
  
  for (int ns = 0; ns<NS; ns++)
  {
	  if (G->Mate[ns] != -1)
	  {
		  Log(" Block %i is identical to %i (reusing T matrices)", ns, G->Mate[ns]);
		  continue;
	  }
	  else
		  Log(" Assembling self contributions to T(%i)...", ns);

	  if (!(G->Surfaces[ns]->IsPEC))
	  {
		  G->RegionMPs[G->Surfaces[ns]->RegionIndices[1]]->UnZero();
		  G->AssembleBEMMatrixBlock(ns, ns, Omega, kBloch, TInt[ns]);
		  G->RegionMPs[G->Surfaces[ns]->RegionIndices[1]]->Zero();
	  };

	  G->RegionMPs[G->Surfaces[ns]->RegionIndices[0]]->UnZero();
	  G->AssembleBEMMatrixBlock(ns, ns, Omega, kBloch, TExt[ns]);
	  G->RegionMPs[G->Surfaces[ns]->RegionIndices[0]]->Zero();
  };
  for (int nr = 0; nr<G->NumRegions; nr++)
	  G->RegionMPs[nr]->UnZero();

  /***************************************************************/
  /* now loop over transformations.                              */
  /* note: 'gtc' stands for 'geometrical transformation complex' */
  /***************************************************************/
  for(int nt=0; nt<SNEQD->NumTransformations; nt++)
   { 
     /*--------------------------------------------------------------*/
     /*- transform the geometry -------------------------------------*/
     /*--------------------------------------------------------------*/
     char *Tag=SNEQD->GTCs[nt]->Tag;
     G->Transform(SNEQD->GTCs[nt]);
     Log(" Computing quantities at geometrical transform %s",Tag);

     /*--------------------------------------------------------------*/
     /* assemble off-diagonal matrix blocks.                         */
     /* note that not all off-diagonal blocks necessarily need to    */
     /* be recomputed for all transformations; this is what the 'if' */
     /* statement here is checking for.                              */
     /*--------------------------------------------------------------*/
     Args->Symmetric=0;
     for(int nb=0, ns=0; ns<NS; ns++)
      for(int nsp=ns+1; nsp<NS; nsp++, nb++)
       if ( nt==0 || G->SurfaceMoved[ns] || G->SurfaceMoved[nsp] )
        G->AssembleBEMMatrixBlock(ns, nsp, Omega, kBloch, U[nb]);
     Log("...SN done with ABMB");

     /*--------------------------------------------------------------*/
     /*- stamp all blocks into the BEM matrix and LU-factorize      -*/
	 /*  in the same loop, store the Exterior blocks into T0         */
     /*--------------------------------------------------------------*/
	 T0->Zero();
     for(int nb=0, ns=0; ns<NS; ns++)
      { 
		 // stamp blocks in the diagonal
        int RowOffset=G->BFIndexOffset[ns];
        M->InsertBlock(TExt[ns], RowOffset, RowOffset);

		// if external region is Exterior, save blocks (F. Ramirez 2019)
		if      (G->Surfaces[ns]->RegionIndices[0] == 0)
			T0->InsertBlock(TExt[ns], RowOffset, RowOffset);
		else if (G->Surfaces[ns]->RegionIndices[1] == 0)
			T0->InsertBlock(TInt[ns], RowOffset, RowOffset);
        
		if( !(G->Surfaces[ns]->IsPEC) )
         M->AddBlock(TInt[ns], RowOffset, RowOffset);

		// stamp blocks in the off diagonal
        for(int nsp=ns+1; nsp<NS; nsp++, nb++)
         { int ColOffset=G->BFIndexOffset[nsp];
           M->InsertBlock(U[nb], RowOffset, ColOffset);
           M->InsertBlockTranspose(U[nb], ColOffset, RowOffset);

		   // if surface region is Exterior, save blocks (F. Ramirez 2019)
		   if (G->Surfaces[nsp]->RegionIndices[0] == 0 || 
			   G->Surfaces[nsp]->RegionIndices[1] == 0){
			   T0->InsertBlock(U[nb], RowOffset, ColOffset);
			   T0->InsertBlockTranspose(U[nb], ColOffset, RowOffset);
		   }
         };
      };
     UndoSCUFFMatrixTransformation(M);
     Log("LU factorizing...");
     M->LUFactorize();
     
     //Log("Computing Green Function gradient from EMT method");
     //gradT0_Matrix(G, SNEQD, Omega, dT0);


     /*--------------------------------------------------------------*/
     /*- compute averge scattering flux                             -*/
     /*--------------------------------------------------------------*/
	 HMatrix *PFTMatrix = SNEQD->PFTMatrix;
	 int NumPFTMethods = SNEQD->NumPFTMethods;

	 AverageScatteringFlux(SNEQD, Omega, PFTMatrix);

	 // Write results to File
	 for (int npm = 0; npm<NumPFTMethods; npm++)
	 {
		 FILE *f = vfopen(SNEQD->SIFluxFileNames[npm], "a");
		 for (int nsd = 0; nsd<NS; nsd++)
		 {
			 fprintf(f, "%s %e ", Tag, real(Omega));
			 if (kBloch) fprintVec(f, kBloch, G->LDim);
			 fprintf(f, "%i%i ", 0, nsd + 1);
			 for (int nq = 0; nq<NUM_AVSCAT; nq++)
				 fprintf(f, "%+.8e ", PFTMatrix->GetEntryD(nsd, nq));
			 fprintf(f, "\n");
		 };
		 fclose(f);

	 };
     /*--------------------------------------------------------------*/
     /* untransform the geometry                                     */
     /*--------------------------------------------------------------*/
     G->UnTransform();
     Log(" ...done!");

  }; // for (nt=0; nt<SNEQD->NumTransformations... )

  /*--------------------------------------------------------------*/
  /*- at the end of the first successful frequency calculation,  -*/
  /*- we dump out the cache to disk, and then tell ourselves not -*/
  /*- to dump the cache to disk again (explain me)               -*/
  /*--------------------------------------------------------------*/
  if ( SNEQD->WriteCache ) 
   { StoreCache( SNEQD->WriteCache );
     SNEQD->WriteCache=0;
   };

}
