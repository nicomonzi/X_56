/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2023
 *
 * Pierangelo Masarati	<pierangelo.masarati@polimi.it>
 * Paolo Mantegazza	<paolo.mantegazza@polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 *
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/* Trave a volumi finiti, con predisposizione per forze di inerzia consistenti
 * e legame cositutivo piezoelettico */


#ifndef BEAM_H
#define BEAM_H

#include <array>

#include "myassert.h"
#include "except.h"

#include "strnode.h"
#include "elem.h"
#include "gravity.h"

#include "constltp.h"

extern const char* psBeamNames[];

/* ... */
class DataManager;
class MBDynParser;

/* Beam - begin */

class Beam
: public InitialAssemblyElem, public GravityOwner {
    friend class AerodynamicBeam;
    friend Elem* ReadBeam(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);
    friend class Beam2;

  public:
    /* Tipi di travi */
    enum Type {
        UNKNOWN = -1,
	ELASTIC = 0,
	VISCOELASTIC,
	PIEZOELECTRICELASTIC,
	PIEZOELECTRICVISCOELASTIC,
	FULLYCOUPLEDPIEZOELECTRICELASTIC,

	LASTBEAMTYPE
    };

	// output mask
	enum {
		OUTPUT_NONE = 0x000U,

		OUTPUT_EP_X = (ToBeOutput::OUTPUT_PRIVATE << 0),
		OUTPUT_EP_R = (ToBeOutput::OUTPUT_PRIVATE << 1),
		OUTPUT_EP_CONFIGURATION = (OUTPUT_EP_X | OUTPUT_EP_R),

		OUTPUT_EP_F = (ToBeOutput::OUTPUT_PRIVATE << 2),
		OUTPUT_EP_M = (ToBeOutput::OUTPUT_PRIVATE << 3),
		OUTPUT_EP_FORCES = (OUTPUT_EP_F | OUTPUT_EP_M ),

		OUTPUT_EP_NU = (ToBeOutput::OUTPUT_PRIVATE << 4),
		OUTPUT_EP_K = (ToBeOutput::OUTPUT_PRIVATE << 5),
		OUTPUT_EP_STRAINS = (OUTPUT_EP_NU | OUTPUT_EP_K),

		OUTPUT_EP_NUP = (ToBeOutput::OUTPUT_PRIVATE << 6),
		OUTPUT_EP_KP = (ToBeOutput::OUTPUT_PRIVATE << 7),
		OUTPUT_EP_STRAIN_RATES = (OUTPUT_EP_NUP | OUTPUT_EP_KP),

		OUTPUT_DEFAULT = (OUTPUT_EP_F | OUTPUT_EP_M),

		OUTPUT_EP_ALL = (ToBeOutput::OUTPUT_PRIVATE_MASK)
	};

protected:
 	static const unsigned int iNumPrivData =
		+3		//  0 ( 1 ->  3) - "e" strain
		+3		//  3 ( 4 ->  6) - "k" curvature
		+3		//  6 ( 7 ->  9) - "F" force
		+3		//  9 (10 -> 12) - "M" moment
		+3		// 12 (13 -> 15) - "X" position
		+3		// 15 (16 -> 18) - "Phi" orientation vector
		+3		// 18 (19 -> 21) - "Omega" angular velocity
		+3		// 21 (22 -> 24) - "eP" strain rate
		+3		// 24 (25 -> 27) - "kP" curvature rate
	;

    static unsigned int iGetPrivDataIdx_int(const char *s,
	ConstLawType::Type type);

	// output flags
	OrientationDescription od;

#ifdef USE_NETCDF
	MBDynNcVar Var_X[2],
		Var_Phi[2],
		Var_F[2],
		Var_M[2],
		Var_Nu[2],
		Var_K[2],
		Var_NuP[2],
		Var_KP[2];
#endif /* USE_NETCDF */

  public:
    class ErrGeneric : public MBDynErrBase {
	public:
		ErrGeneric(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};

  public:
    enum Section { S_I = 0, SII = 1, NUMSEZ = 2 };
    enum NodeName { NODE1 = 0, NODE2 = 1, NODE3 = 2, NUMNODES = 3 };
    enum Deformations { STRAIN = 0, CURVAT = 1, NUMDEFORM = 2 };

  protected:
    /* Puntatori ai nodi */
    std::array<const StructNode*, NUMNODES> pNode;

    /* Offset dei nodi */
    const std::array<Vec3, NUMNODES> f;
    std::array<Vec3, NUMNODES> fRef;
    const std::array<Mat3x3, NUMNODES> RNode;

    /* Matrice di rotazione delle sezioni - non sono const perche' vengono
     * aggiornate ad ogni iterazione */
    std::array<Mat3x3, NUMSEZ> R;
    std::array<Mat3x3, NUMSEZ> RRef;
    std::array<Mat3x3, NUMSEZ> RPrev;

    /* Constitutive laws*/
    std::array<ConstitutiveLaw6D*, NUMSEZ> pD;

    /* Reference constitutive laws */
    std::array<Mat6x6, NUMSEZ> DRef;

    /* Per forze d'inerzia consistenti: */
    const bool bConsistentInertia;

    const doublereal dMass_I;
    const Vec3 S0_I;
    const Mat3x3 J0_I;

    const doublereal dMassII;
    const Vec3 S0II;
    const Mat3x3 J0II;

    /* Velocita' angolare delle sezioni */
    std::array<Vec3, NUMSEZ> Omega;
    std::array<Vec3, NUMSEZ> OmegaRef;

    /* Dati temporanei che vengono passati da AssRes ad AssJac */
    std::array<Vec6, NUMSEZ> Az;
    std::array<Vec6, NUMSEZ> AzRef;
    std::array<Vec6, NUMSEZ> AzLoc;
    std::array<Vec6, NUMSEZ> DefLoc;
    std::array<Vec6, NUMSEZ> DefLocRef;
    std::array<Vec6, NUMSEZ> DefLocPrev;

    // NOTE: Moved to Beam from ViscoElasticBeam for output purposes
    std::array<Vec6, NUMSEZ> DefPrimeLoc;

    std::array<Vec3, NUMSEZ> p;
    std::array<Vec3, NUMSEZ> g;
    std::array<Vec3, NUMSEZ> L0;
    std::array<Vec3, NUMSEZ> L;

    std::array<Vec3, NUMSEZ> LRef;

    std::array<doublereal, NUMSEZ> dsdxi;

    /* Is first res? */
    bool bFirstRes;

    /* Funzioni di servizio */
    static Vec3
    InterpState(const Vec3& v1,
                const Vec3& v2,
		const Vec3& v3,
		enum Section Sec);
    Vec3
    InterpDeriv(const Vec3& v1,
                const Vec3& v2,
		const Vec3& v3,
		enum Section Sec);

    /* Funzioni di calcolo delle matrici */
    virtual void
    AssStiffnessMat(FullSubMatrixHandler& WMA,
                    FullSubMatrixHandler& WMB,
		    doublereal dCoef,
		    const VectorHandler& XCurr,
		    const VectorHandler& XPrimeCurr);

    /* Per le beam che aggiungono qualcosa alle az. interne */
    virtual void
    AddInternalForces(Vec6& /* AzLoc */ , unsigned int /* iSez */ ) {
        NO_OP;
    };
     
    virtual void
    AssStiffnessVec(SubVectorHandler& WorkVec,
                    doublereal dCoef,
		    const VectorHandler& XCurr,
		    const VectorHandler& XPrimeCurr);
     
    virtual void
    AssInertiaMat(FullSubMatrixHandler& /* WMA */ ,
                  FullSubMatrixHandler& /* WMB */ ,
		  doublereal /* dCoef */ ,
		  const VectorHandler& /* XCurr */ ,
		  const VectorHandler& /* XPrimeCurr */ ) {
        NO_OP;
    };

    virtual void
    AssInertiaVec(SubVectorHandler& /* WorkVec */ ,
                  doublereal /* dCoef */ ,
		  const VectorHandler& /* XCurr */ ,
		  const VectorHandler& /* XPrimeCurr */ ) {
        NO_OP;
    };

    /* Inizializza le derivate delle funzioni di forma
     * e calcola le deformazioni geometriche iniziali */
    virtual void DsDxi(void);

    /* Calcola la velocita' angolare delle sezioni a partire da quelle
     * dei nodi; per ora lo fa in modo brutale, ma si puo' fare anche
     * in modo consistente */
    virtual void Omega0(void);

    /* Funzione interna di restart */
    virtual std::ostream& Restart_(std::ostream& out) const;

    /* Inizializza i dati */
    void Init(void);

  public:
    /* Costruttore normale */
    Beam(unsigned int uL,
	 const StructNode* pN1, const StructNode* pN2, const StructNode* pN3,
	 const Vec3& F1, const Vec3& F2, const Vec3& F3,
	 const Mat3x3& R1, const Mat3x3& R2, const Mat3x3& R3,
	 const Mat3x3& r_I, const Mat3x3& rII,
	 ConstitutiveLaw6D* const pD_I, ConstitutiveLaw6D* const pDII,
	 OrientationDescription ood,
	 flag fOut);

    /* Costruttore per la trave con forze d'inerzia consistenti */
    Beam(unsigned int uL,
	 const StructNode* pN1, const StructNode* pN2, const StructNode* pN3,
	 const Vec3& F1, const Vec3& F2, const Vec3& F3,
	 const Mat3x3& R1, const Mat3x3& R2, const Mat3x3& R3,
	 const Mat3x3& r_I, const Mat3x3& rII,
	 ConstitutiveLaw6D* const pD_I, ConstitutiveLaw6D* const pDII,
	 doublereal dM_I,
	 const Vec3& s0_I, const Mat3x3& j0_I,
	 doublereal dMII,
	 const Vec3& s0II, const Mat3x3& j0II,
	 OrientationDescription ood,
	 flag fOut);

    /* Distruttore banale */
    virtual ~Beam(void);

    /* Tipo di trave */
    virtual Beam::Type GetBeamType(void) const {
        return Beam::ELASTIC;
    };

    /* Tipo di elemento */
    virtual Elem::Type GetElemType(void) const override {
        return Elem::BEAM;
    };
    
    /* Deformable element */
    virtual bool bIsDeformable() const override {
        return true;
    };

    /* Contributo al file di restart */
    virtual std::ostream& Restart(std::ostream& out) const override;

    virtual void Restart(RestartData& oData, RestartData::RestartAction eAction) override;     
     
    virtual void
    AfterConvergence(const VectorHandler& X, const VectorHandler& XP) override;

    /* Inverse Dynamics */
    virtual void
    AfterConvergence(const VectorHandler& X, const VectorHandler& XP,
    		const VectorHandler& XPP) override;

    /* funzioni proprie */

    /* Dimensioni del workspace; sono 36 righe perche' se genera anche le
     * forze d'inerzia consistenti deve avere accesso alle righe di definizione
     * della quantita' di moto */
    virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const override {
        if (bConsistentInertia) {
	    *piNumRows = 36;
        } else {
	    *piNumRows = 18;
        }

        *piNumCols = 18;
    };

    /* Settings iniziali, prima della prima soluzione */
    void SetValue(DataManager *pDM,
		    VectorHandler& /* X */ , VectorHandler& /* XP */ ,
		    SimulationEntity::Hints *ph = 0) override;

    /* Prepara i parametri di riferimento dopo la predizione */
    virtual void
    AfterPredict(VectorHandler& /* X */ , VectorHandler& /* XP */ ) override;

    /* assemblaggio residuo */
    virtual SubVectorHandler&
    AssRes(SubVectorHandler& WorkVec,
           doublereal dCoef,
	   const VectorHandler& XCurr,
	   const VectorHandler& XPrimeCurr) override;

    /* Inverse Dynamics: */
    virtual SubVectorHandler&
    AssRes(SubVectorHandler& WorkVec,
	   const VectorHandler&  XCurr ,
	   const VectorHandler&  XPrimeCurr ,
	   const VectorHandler&  XPrimePrimeCurr ,
	   InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS) override;

    /* inverse dynamics capable element */
    virtual bool bInverseDynamics(void) const override;

    /* assemblaggio jacobiano */
    virtual VariableSubMatrixHandler&
    AssJac(VariableSubMatrixHandler& WorkMat,
	   doublereal dCoef,
	   const VectorHandler& XCurr,
	   const VectorHandler& XPrimeCurr) override;

    /* assemblaggio matrici per autovalori */
    void
    AssMats(VariableSubMatrixHandler& WorkMatA,
	   VariableSubMatrixHandler& WorkMatB,
	   const VectorHandler& XCurr,
	   const VectorHandler& XPrimeCurr) override;

    virtual void OutputPrepare(OutputHandler &OH) override;

    /* output; si assume che ogni tipo di elemento sappia, attraverso
     * l'OutputHandler, dove scrivere il proprio output */
    virtual void Output(OutputHandler& OH) const override;

#if 0
    /* Output di un modello NASTRAN equivalente nella configurazione corrente */
    virtual void Output_pch(std::ostream& out) const;
#endif

    /* Funzioni proprie tipiche dei vincoli, usate durante l'assemblaggio
     * iniziale. Le travi non sono vincoli (o meglio, non vengono considerate
     * tali), ma richiedono comunque assemblaggio iniziale. Innanzitutto
     * le coordinate dei nodi che sono fornite vengono usate per calcolare
     * la lunghezza della trave e per la metrica della trave stessa.
     * Inoltre durante l'assemblaggio la trave esplica la propria
     * deformabilita' quando i vincoli tentano di muovere i nodi. */

    /* Numero di gradi di liberta' definiti durante l'assemblaggio iniziale
     * e' dato dai gradi di liberta' soliti piu' le loro derivate necessarie;
     * tipicamente per un vincolo di posizione il numero di dof raddoppia, in
     * quanto vengono derivate tutte le equazioni, mentre per un vincolo di
     * velocita' rimane inalterato. Sono possibili casi intermedi per vincoli
     * misti di posizione e velocita' */
    virtual unsigned int iGetInitialNumDof(void) const override {
        return 0;
    };

    /* Dimensione del workspace durante l'assemblaggio iniziale. Occorre tener
     * conto del numero di dof che l'elemento definisce in questa fase e dei
     * dof dei nodi che vengono utilizzati. Sono considerati dof indipendenti
     * la posizione e la velocita' dei nodi */
    virtual void
    InitialWorkSpaceDim(integer* piNumRows,
			integer* piNumCols) const override {
        *piNumRows = 18;
        *piNumCols = 18;
    };

    /* Setta il valore iniziale delle proprie variabili */
    virtual void SetInitialValue(VectorHandler& /* X */ ) {
        NO_OP;
    };

    /* Contributo allo jacobiano durante l'assemblaggio iniziale */
    virtual VariableSubMatrixHandler&
    InitialAssJac(VariableSubMatrixHandler& WorkMat,
                  const VectorHandler& XCurr) override;

    /* Contributo al residuo durante l'assemblaggio iniziale */
    virtual SubVectorHandler&
    InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr) override;

    /* Accesso ai dati privati */
    virtual unsigned int iGetNumPrivData(void) const override;
    virtual unsigned int iGetPrivDataIdx(const char *s) const override;
    virtual doublereal dGetPrivData(unsigned int i) const override;

    /* Accesso ai nodi */
    virtual const StructNode* pGetNode(unsigned int i) const;

    /******** PER IL SOLUTORE PARALLELO *********/
    /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
     * utile per l'assemblaggio della matrice di connessione fra i dofs */
    virtual void
    GetConnectedNodes(std::vector<const Node *>& connectedNodes) const override {
        connectedNodes.resize(NUMNODES);
        for (int i = 0; i < NUMNODES; i++) {
            connectedNodes[i] = pNode[i];
        }
    };
    /**************************************************/

};

/* Beam - end */


/* ViscoElasticBeam - begin */

class ViscoElasticBeam : public Beam {
  protected:

    /* Derivate di deformazioni e curvature */
     std::array<Vec3, NUMSEZ> LPrime;
     std::array<Vec3, NUMSEZ> gPrime;

     std::array<Vec3, NUMSEZ> LPrimeRef;

    // NOTE: Moved to Beam from ViscoElasticBeam for output purposes
    // Vec6 DefPrimeLoc[NUMSEZ];

     std::array<Vec6, NUMSEZ> DefPrimeLocRef;

     std::array<Mat6x6, NUMSEZ> ERef;

    /* Funzioni di calcolo delle matrici */
    virtual void
    AssStiffnessMat(FullSubMatrixHandler& WMA,
                    FullSubMatrixHandler& WMB,
		    doublereal dCoef,
		    const VectorHandler& XCurr,
		    const VectorHandler& XPrimeCurr) override;

    virtual void
    AssStiffnessVec(SubVectorHandler& WorkVec,
                    doublereal dCoef,
		    const VectorHandler& XCurr,
		    const VectorHandler& XPrimeCurr) override;

    /* Inizializza i dati */
    void Init(void);

  public:
    /* Costruttore normale */
    ViscoElasticBeam(unsigned int uL,
	             const StructNode* pN1,
		     const StructNode* pN2,
		     const StructNode* pN3,
	             const Vec3& F1,
		     const Vec3& F2,
		     const Vec3& F3,
		     const Mat3x3& R1,
		     const Mat3x3& R2,
		     const Mat3x3& R3,
	             const Mat3x3& r_I,
		     const Mat3x3& rII,
	             ConstitutiveLaw6D* const pD_I,
		     ConstitutiveLaw6D* const pDII,
		     OrientationDescription ood,
		     flag fOut);

    /* Costruttore per la trave con forze d'inerzia consistenti */
    ViscoElasticBeam(unsigned int uL,
                     const StructNode* pN1,
		     const StructNode* pN2,
		     const StructNode* pN3,
		     const Vec3& F1,
		     const Vec3& F2,
		     const Vec3& F3,
		     const Mat3x3& R1,
		     const Mat3x3& R2,
		     const Mat3x3& R3,
		     const Mat3x3& r_I,
		     const Mat3x3& rII,
		     ConstitutiveLaw6D* const pD_I,
		     ConstitutiveLaw6D* const pDII,
		     doublereal dM_I,
		     const Vec3& s0_I,
		     const Mat3x3& j0_I,
		     doublereal dMII,
		     const Vec3& s0II,
		     const Mat3x3& j0II,
		     OrientationDescription ood,
		     flag fOut);

    /* Distruttore banale */
    virtual ~ViscoElasticBeam(void) {
        NO_OP;
    };

    /* Tipo di trave */
    virtual Beam::Type GetBeamType(void) const override {
        return Beam::VISCOELASTIC;
    };

    /* Settings iniziali, prima della prima soluzione */
    void SetValue(DataManager *pDM,
                  VectorHandler& /* X */ , VectorHandler& /* XP */ ,
                  SimulationEntity::Hints *ph = 0) override;

    /* Prepara i parametri di riferimento dopo la predizione */
    virtual void
    AfterPredict(VectorHandler& /* X */ , VectorHandler& /* XP */ ) override;

    virtual void
    AfterConvergence(const VectorHandler& X, const VectorHandler& XP) override;

    /* Inverse Dynamics */
    virtual void
    AfterConvergence(const VectorHandler& X, const VectorHandler& XP,
    		const VectorHandler& XPP) override;

    virtual doublereal dGetPrivData(unsigned int i) const override;

    virtual void Restart(RestartData& oData, RestartData::RestartAction eAction) override;
};

/* ViscoElasticBeam - end */

extern void
ReadBeamCustomOutput(DataManager* pDM, MBDynParser& HP, unsigned int uLabel,
	Beam::Type BT, unsigned& uFlags, OrientationDescription& od);

extern void
ReadOptionalBeamCustomOutput(DataManager* pDM, MBDynParser& HP, unsigned int uLabel,
	Beam::Type BT, unsigned& uFlags, OrientationDescription& od);

extern Elem*
ReadBeam(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);

#endif /* BEAM_H */
