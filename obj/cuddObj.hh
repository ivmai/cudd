/**CHeaderFile***************************************************************

  FileName    [cuddObj.hh]

  PackageName [cudd]

  Synopsis    [Class definitions for C++ object-oriented encapsulation of
  CUDD.]

  Description [Class definitions for C++ object-oriented encapsulation of
  CUDD.]

  SeeAlso     []

  Author      [Fabio Somenzi]

  Copyright [This file was created at the University of Colorado at
  Boulder.  The University of Colorado at Boulder makes no warranty
  about the suitability of this software for any purpose.  It is
  presented on an AS IS basis.]

  Revision    [$Id: cuddObj.hh,v 1.5 2004/01/01 07:04:22 fabio Exp fabio $]

******************************************************************************/

#ifndef _CPPCUDD
#define _CPPCUDD

/*---------------------------------------------------------------------------*/
/* Nested includes                                                           */
/*---------------------------------------------------------------------------*/

#include <string>
#include <iostream>
#include "util.h"
#include "cudd.h"

using std::cout;
using std::cerr;
using std::endl;
using std::hex;
using std::string;

/*---------------------------------------------------------------------------*/
/* Type definitions                                                          */
/*---------------------------------------------------------------------------*/

class ADD;
class BDD;
class ZDD;
class Cudd;
class BDDvector;
class ADDvector;
class ZDDvector;

typedef void (*PFC)(string);	// error function type

/*---------------------------------------------------------------------------*/
/* Class definitions                                                         */
/*---------------------------------------------------------------------------*/


/**Class***********************************************************************

  Synopsis     [Base class for all decision diagrams in CUDD.]

  Description  []

  SeeAlso      [Cudd ABDD ADD BDD ZDD]

******************************************************************************/
class DD {
    friend class ABDD;
    friend class ADD;
    friend class BDD;
    friend class ZDD;
    Cudd *ddMgr;
    DdNode *node;
    inline DdManager * checkSameManager(const DD &other);
    inline void checkReturnValue(const DdNode *result);
    inline void checkReturnValue(const int result, const int expected = 1);
public:
    DD(Cudd *ddManager, DdNode *ddNode);
    DD();
    DD(const DD &from);
    Cudd *manager();
    inline DdNode * getNode();
    int nodeCount();
    unsigned int NodeReadIndex();

}; // DD


/**Class***********************************************************************

  Synopsis     [Class for ADDs and BDDs.]

  Description  []

  SeeAlso      [Cudd ADD BDD]

******************************************************************************/
class ABDD : public DD {
    friend class BDD;
    friend class ADD;
    friend class Cudd;
public:
    ABDD(Cudd *bddManager, DdNode *bddNode);
    ABDD();
    ABDD(const ABDD &from);
    virtual ~ABDD();
    int operator==(const ABDD &other);
    int operator!=(const ABDD &other);
    void print(int nvars, int verbosity = 1);
    DdApaNumber ApaCountMinterm(int nvars, int * digits);
    void ApaPrintMinterm(int nvars, FILE * fp = stdout);
    void EpdPrintMinterm(int nvars, FILE * fp = stdout);
    BDD FindEssential();
    void PrintTwoLiteralClauses(char ** names, FILE * fp = stdout);
    BDD ShortestPath(int * weight, int * support, int * length);
    BDD LargestCube(int * length);
    int ShortestLength(int * weight);
    int EquivDC(ABDD G, ABDD D);
    double * CofMinterm();
    void PrintMinterm();
    double CountMinterm(int nvars);
    double CountPath();
    BDD Support();
    int SupportSize();
    void ClassifySupport(ABDD g, BDD* common, BDD* onlyF, BDD* onlyG);
    int CountLeaves();
    DdGen * FirstCube(int ** cube, CUDD_VALUE_TYPE * value);
    double Density(int nvars);

}; // ABDD


/**Class***********************************************************************

  Synopsis     [Class for BDDs.]

  Description  []

  SeeAlso      [Cudd]

******************************************************************************/
class BDD : public ABDD {
    friend class Cudd;
public:
    BDD(Cudd *bddManager, DdNode *bddNode);
    BDD();
    BDD(const BDD &from);
    int operator==(const BDD& other);
    int operator!=(const BDD& other);
    BDD operator=(const BDD& right);
    int operator<=(const BDD& other);
    int operator>=(const BDD& other);
    int operator<(const BDD& other);
    int operator>(const BDD& other);
    BDD operator!();
    BDD operator~();
    BDD operator*(const BDD& other);
    BDD operator*=(const BDD& other);
    BDD operator&(const BDD& other);
    BDD operator&=(const BDD& other);
    BDD operator+(const BDD& other);
    BDD operator+=(const BDD& other);
    BDD operator|(const BDD& other);
    BDD operator|=(const BDD& other);
    BDD operator^(const BDD& other);
    BDD operator^=(const BDD& other);
    BDD operator-(const BDD& other);
    BDD operator-=(const BDD& other);
    BDD AndAbstract(BDD g, BDD cube);
    BDD AndAbstractLimit(BDD g, BDD cube, unsigned int limit);
    BDD UnderApprox(
      int numVars,
      int threshold = 0,
      int safe = 0,
      double quality = 1.0);
    BDD OverApprox(
      int numVars,
      int threshold = 0,
      int safe = 0,
      double quality = 1.0);
    BDD RemapUnderApprox(int numVars, int threshold = 0, double quality = 1.0);
    BDD RemapOverApprox(int numVars, int threshold = 0, double quality = 1.0);
    BDD ExistAbstract(BDD cube);
    BDD XorExistAbstract(BDD g, BDD cube);
    BDD UnivAbstract(BDD cube);
    BDD BooleanDiff(int x);
    int VarIsDependent(BDD var);
    double Correlation(BDD g);
    double CorrelationWeights(BDD g, double * prob);
    BDD Ite(BDD g, BDD h);
    BDD IteConstant(BDD g, BDD h);
    BDD Intersect(BDD g);
    BDD And(BDD g);
    BDD AndLimit(BDD g, unsigned int limit);
    BDD Or(BDD g);
    BDD Nand(BDD g);
    BDD Nor(BDD g);
    BDD Xor(BDD g);
    BDD Xnor(BDD g);
    int Leq(BDD g);
    ADD Add();
    BDD Transfer(Cudd& destination);
    BDD ClippingAnd(BDD g, int maxDepth, int direction);
    BDD ClippingAndAbstract(BDD g, BDD cube, int maxDepth, int direction);
    BDD Cofactor(BDD g);
    BDD Compose(BDD g, int v);
    BDD Permute(int * permut);
    BDD SwapVariables(BDDvector x, BDDvector y);
    BDD AdjPermuteX(BDDvector x);
    BDD VectorCompose(BDDvector vector);
    void ApproxConjDecomp(BDD* g, BDD* h);
    void ApproxDisjDecomp(BDD* g, BDD* h);
    void IterConjDecomp(BDD* g, BDD* h);
    void IterDisjDecomp(BDD* g, BDD* h);
    void GenConjDecomp(BDD* g, BDD* h);
    void GenDisjDecomp(BDD* g, BDD* h);
    void VarConjDecomp(BDD* g, BDD* h);
    void VarDisjDecomp(BDD* g, BDD* h);
    int IsVarEssential(int id, int phase);
    BDD Constrain(BDD c);
    BDD Restrict(BDD c);
    BDD NPAnd(BDD g);
    BDDvector ConstrainDecomp();
    BDDvector CharToVect();
    BDD LICompaction(BDD c);
    BDD Squeeze(BDD u);
    BDD Minimize(BDD c);
    BDD SubsetCompress(int nvars, int threshold);
    BDD SupersetCompress(int nvars, int threshold);
    BDD LiteralSetIntersection(BDD g);
    BDD PrioritySelect(BDDvector x, BDDvector y, BDDvector z, BDD Pi,
		       DD_PRFP Pifunc);
    BDD CProjection(BDD Y);
    int MinHammingDist(int *minterm, int upperBound);
    BDD Eval(int * inputs);
    BDD Decreasing(int i);
    BDD Increasing(int i);
    BDD SolveEqn(BDD Y, BDD* G, int ** yIndex, int n);
    BDD VerifySol(BDD* G, int * yIndex, int n);
    BDD SplitSet(BDDvector xVars, double m);
    BDD SubsetHeavyBranch(int numVars, int threshold);
    BDD SupersetHeavyBranch(int numVars, int threshold);
    BDD SubsetShortPaths(int numVars, int threshold, int hardlimit);
    BDD SupersetShortPaths(int numVars, int threshold, int hardlimit);
    void PrintCover();
    void PrintCover(const BDD& u);
    int EstimateCofactor(int i, int phase);
    int EstimateCofactorSimple(int i);
    void PickOneCube(char * string);
    BDD PickOneMinterm(BDDvector vars);
    DdGen * FirstNode(BDD* fnode);
    BDD zddIsop(BDD U, ZDD* zdd_I);
    BDD Isop(BDD U);
    ZDD PortToZdd();

}; // BDD


/**Class***********************************************************************

  Synopsis     [Class for ADDs.]

  Description  []

  SeeAlso      [Cudd]

******************************************************************************/
class ADD : public ABDD {
    friend class Cudd;
public:
    ADD(Cudd *bddManager, DdNode *bddNode);
    ADD();
    ADD(const ADD &from);
    int operator==(const ADD& other);
    int operator!=(const ADD& other);
    ADD operator=(const ADD& right);
    // Relational operators
    int operator<=(const ADD& other);
    int operator>=(const ADD& other);
    int operator<(const ADD& other);
    int operator>(const ADD& other);
    // Arithmetic operators
    ADD operator-();
    ADD operator*(const ADD& other);
    ADD operator*=(const ADD& other);
    ADD operator+(const ADD& other);
    ADD operator+=(const ADD& other);
    ADD operator-(const ADD& other);
    ADD operator-=(const ADD& other);
    // Logical operators
    ADD operator~();
    ADD operator&(const ADD& other);
    ADD operator&=(const ADD& other);
    ADD operator|(const ADD& other);
    ADD operator|=(const ADD& other);
    ADD ExistAbstract(ADD cube);
    ADD UnivAbstract(ADD cube);
    ADD OrAbstract(ADD cube);
    ADD Plus(ADD g);
    ADD Times(ADD g);
    ADD Threshold(ADD g);
    ADD SetNZ(ADD g);
    ADD Divide(ADD g);
    ADD Minus(ADD g);
    ADD Minimum(ADD g);
    ADD Maximum(ADD g);
    ADD OneZeroMaximum(ADD g);
    ADD Diff(ADD g);
    ADD Agreement(ADD g);
    ADD Or(ADD g);
    ADD Nand(ADD g);
    ADD Nor(ADD g);
    ADD Xor(ADD g);
    ADD Xnor(ADD g);
    ADD Log();
    ADD FindMax();
    ADD FindMin();
    ADD IthBit(int bit);
    ADD ScalarInverse(ADD epsilon);
    ADD Ite(ADD g, ADD h);
    ADD IteConstant(ADD g, ADD h);
    ADD EvalConst(ADD g);
    int Leq(ADD g);
    ADD Cmpl();
    ADD Negate();
    ADD RoundOff(int N);
    BDD BddThreshold(CUDD_VALUE_TYPE value);
    BDD BddStrictThreshold(CUDD_VALUE_TYPE value);
    BDD BddInterval(CUDD_VALUE_TYPE lower, CUDD_VALUE_TYPE upper);
    BDD BddIthBit(int bit);
    BDD BddPattern();
    ADD Cofactor(ADD g);
    ADD Compose(ADD g, int v);
    ADD Permute(int * permut);
    ADD SwapVariables(ADDvector x, ADDvector y);
    ADD VectorCompose(ADDvector vector);
    ADD NonSimCompose(ADDvector vector);
    ADD Constrain(ADD c);
    ADD Restrict(ADD c);
    ADD MatrixMultiply(ADD B, ADDvector z);
    ADD TimesPlus(ADD B, ADDvector z);
    ADD Triangle(ADD g, ADDvector z);
    ADD Eval(int * inputs);
    int EqualSupNorm(ADD g, CUDD_VALUE_TYPE tolerance, int pr);

}; // ADD


/**Class***********************************************************************

  Synopsis     [Class for ZDDs.]

  Description  []

  SeeAlso      [Cudd]

******************************************************************************/
class ZDD : public DD {
    friend class Cudd;
public:
    ZDD(Cudd *bddManager, DdNode *bddNode);
    ZDD();
    ZDD(const ZDD &from);
    ~ZDD();
    ZDD operator=(const ZDD& right);
    int operator==(const ZDD& other);
    int operator!=(const ZDD& other);
    int operator<=(const ZDD& other);
    int operator>=(const ZDD& other);
    int operator<(const ZDD& other);
    int operator>(const ZDD& other);
    void print(int nvars, int verbosity = 1);
    ZDD operator*(const ZDD& other);
    ZDD operator*=(const ZDD& other);
    ZDD operator&(const ZDD& other);
    ZDD operator&=(const ZDD& other);
    ZDD operator+(const ZDD& other);
    ZDD operator+=(const ZDD& other);
    ZDD operator|(const ZDD& other);
    ZDD operator|=(const ZDD& other);
    ZDD operator-(const ZDD& other);
    ZDD operator-=(const ZDD& other);
    int Count();
    double CountDouble();
    ZDD Product(ZDD g);
    ZDD UnateProduct(ZDD g);
    ZDD WeakDiv(ZDD g);
    ZDD Divide(ZDD g);
    ZDD WeakDivF(ZDD g);
    ZDD DivideF(ZDD g);
    double CountMinterm(int path);
    BDD PortToBdd();
    ZDD Ite(ZDD g, ZDD h);
    ZDD Union(ZDD Q);
    ZDD Intersect(ZDD Q);
    ZDD Diff(ZDD Q);
    ZDD DiffConst(ZDD Q);
    ZDD Subset1(int var);
    ZDD Subset0(int var);
    ZDD Change(int var);
    void PrintMinterm();
    void PrintCover();

}; // ZDD


/**Class***********************************************************************

  Synopsis     [Class for CUDD managers.]

  Description  []

  SeeAlso      [DD]

******************************************************************************/
class Cudd {
    friend class DD;
    friend class ABDD;
    friend class ADD;
    friend class BDD;
    friend class ZDD;
    struct capsule {
	DdManager *manager;
	PFC errorHandler;
	int verbose;
	int ref;
    };
    capsule *p;
public:
    Cudd(
      unsigned int numVars = 0,
      unsigned int numVarsZ = 0,
      unsigned int numSlots = CUDD_UNIQUE_SLOTS,
      unsigned int cacheSize = CUDD_CACHE_SLOTS,
      unsigned long maxMemory = 0);
    Cudd(Cudd& x);
    ~Cudd();
    PFC setHandler(PFC newHandler);
    PFC getHandler();
    DdManager *getManager() {return p->manager;}
    inline void makeVerbose() {p->verbose = 1;}
    inline void makeTerse() {p->verbose = 0;}
    inline int isVerbose() {return p->verbose;}
    inline void checkReturnValue(const DdNode *result);
    inline void checkReturnValue(const int result);
    Cudd& operator=(const Cudd& right);
    void info();
    BDD bddVar();
    BDD bddVar(int index);
    BDD bddOne();
    BDD bddZero();
    ADD addVar();
    ADD addVar(int index);
    ADD addOne();
    ADD addZero();
    ADD constant(CUDD_VALUE_TYPE c);
    ADD plusInfinity();
    ADD minusInfinity();
    ZDD zddVar(int index);
    ZDD zddOne(int i);
    ZDD zddZero();
    ADD addNewVarAtLevel(int level);
    BDD bddNewVarAtLevel(int level);
    void zddVarsFromBddVars(int multiplicity);
    void AutodynEnable(Cudd_ReorderingType method);
    void AutodynDisable();
    int ReorderingStatus(Cudd_ReorderingType * method);
    void AutodynEnableZdd(Cudd_ReorderingType method);
    void AutodynDisableZdd();
    int ReorderingStatusZdd(Cudd_ReorderingType * method);
    int zddRealignmentEnabled();
    void zddRealignEnable();
    void zddRealignDisable();
    int bddRealignmentEnabled();
    void bddRealignEnable();
    void bddRealignDisable();
    ADD background();
    void SetBackground(ADD bg);
    unsigned int ReadCacheSlots();
    double ReadCacheUsedSlots();
    double ReadCacheLookUps();
    double ReadCacheHits();
    unsigned int ReadMinHit();
    void SetMinHit(unsigned int hr);
    unsigned int ReadLooseUpTo();
    void SetLooseUpTo(unsigned int lut);
    unsigned int ReadMaxCache();
    unsigned int ReadMaxCacheHard();
    void SetMaxCacheHard(unsigned int mc);
    int ReadSize();
    int ReadZddSize();
    unsigned int ReadSlots();
    unsigned int ReadKeys();
    unsigned int ReadDead();
    unsigned int ReadMinDead();
    int ReadReorderings();
    long ReadReorderingTime();
    int ReadGarbageCollections();
    long ReadGarbageCollectionTime();
    int ReadSiftMaxVar();
    void SetSiftMaxVar(int smv);
    int ReadSiftMaxSwap();
    void SetSiftMaxSwap(int sms);
    double ReadMaxGrowth();
    void SetMaxGrowth(double mg);
    MtrNode * ReadTree();
    void SetTree(MtrNode * tree);
    void FreeTree();
    MtrNode * ReadZddTree();
    void SetZddTree(MtrNode * tree);
    void FreeZddTree();
    int ReadPerm(int i);
    int ReadPermZdd(int i);
    int ReadInvPerm(int i);
    int ReadInvPermZdd(int i);
    BDD ReadVars(int i);
    CUDD_VALUE_TYPE ReadEpsilon();
    void SetEpsilon(CUDD_VALUE_TYPE ep);
    Cudd_AggregationType ReadGroupcheck();
    void SetGroupcheck(Cudd_AggregationType gc);
    int GarbageCollectionEnabled();
    void EnableGarbageCollection();
    void DisableGarbageCollection();
    int DeadAreCounted();
    void TurnOnCountDead();
    void TurnOffCountDead();
    int ReadRecomb();
    void SetRecomb(int recomb);
    int ReadSymmviolation();
    void SetSymmviolation(int symmviolation);
    int ReadArcviolation();
    void SetArcviolation(int arcviolation);
    int ReadPopulationSize();
    void SetPopulationSize(int populationSize);
    int ReadNumberXovers();
    void SetNumberXovers(int numberXovers);
    unsigned long ReadMemoryInUse();
    long ReadPeakNodeCount();
    long ReadNodeCount();
    long zddReadNodeCount();
    void AddHook(DD_HFP f, Cudd_HookType where);
    void RemoveHook(DD_HFP f, Cudd_HookType where);
    int IsInHook(DD_HFP f, Cudd_HookType where);
    void EnableReorderingReporting();
    void DisableReorderingReporting();
    int ReorderingReporting();
    int ReadErrorCode();
    void ClearErrorCode();
    FILE *ReadStdout();
    void SetStdout(FILE *);
    FILE *ReadStderr();
    void SetStderr(FILE *);
    unsigned int ReadNextReordering();
    double ReadSwapSteps();
    unsigned int ReadMaxLive();
    void SetMaxLive(unsigned int);
    unsigned long ReadMaxMemory();
    void SetMaxMemory(unsigned long);
    int bddBindVar(int);
    int bddUnbindVar(int);
    int bddVarIsBound(int);
    ADD Walsh(ADDvector x, ADDvector y);
    ADD addResidue(int n, int m, int options, int top);
    int ApaNumberOfDigits(int binaryDigits);
    DdApaNumber NewApaNumber(int digits);
    void ApaCopy(int digits, DdApaNumber source, DdApaNumber dest);
    DdApaDigit ApaAdd(int digits, DdApaNumber a, DdApaNumber b, DdApaNumber sum);
    DdApaDigit ApaSubtract(int digits, DdApaNumber a, DdApaNumber b, DdApaNumber diff);
    DdApaDigit ApaShortDivision(int digits, DdApaNumber dividend, DdApaDigit divisor, DdApaNumber quotient);
    void ApaShiftRight(int digits, DdApaDigit in, DdApaNumber a, DdApaNumber b);
    void ApaSetToLiteral(int digits, DdApaNumber number, DdApaDigit literal);
    void ApaPowerOfTwo(int digits, DdApaNumber number, int power);
    void ApaPrintHex(FILE * fp, int digits, DdApaNumber number);
    void ApaPrintDecimal(FILE * fp, int digits, DdApaNumber number);
    void DebugCheck();
    void CheckKeys();
    MtrNode * MakeTreeNode(unsigned int low, unsigned int size, unsigned int type);
    // void Harwell(FILE * fp, ADD* E, ADD** x, ADD** y, ADD** xn, ADD** yn_, int * nx, int * ny, int * m, int * n, int bx, int sx, int by, int sy, int pr);
    void PrintLinear();
    int ReadLinear(int x, int y);
    BDD Xgty(BDDvector z, BDDvector x, BDDvector y);
    BDD Xeqy(BDDvector x, BDDvector y);
    ADD Xeqy(ADDvector x, ADDvector y);
    BDD Dxygtdxz(BDDvector x, BDDvector y, BDDvector z);
    BDD Dxygtdyz(BDDvector x, BDDvector y, BDDvector z);
    ADD Hamming(ADDvector xVars, ADDvector yVars);
    // void Read(FILE * fp, ADD* E, ADD** x, ADD** y, ADD** xn, ADD** yn_, int * nx, int * ny, int * m, int * n, int bx, int sx, int by, int sy);
    // void Read(FILE * fp, BDD* E, BDD** x, BDD** y, int * nx, int * ny, int * m, int * n, int bx, int sx, int by, int sy);
    void ReduceHeap(Cudd_ReorderingType heuristic, int minsize);
    void ShuffleHeap(int * permutation);
    void SymmProfile(int lower, int upper);
    unsigned int Prime(unsigned int pr);
    int SharingSize(DD* nodes, int n);
    BDD bddComputeCube(BDD * vars, int * phase, int n);
    ADD addComputeCube(ADD * vars, int * phase, int n);
    int NextNode(DdGen * gen, BDD * nnode);
    BDD IndicesToCube(int * array, int n);
    void PrintVersion(FILE * fp);
    double AverageDistance();
    long Random();
    void Srandom(long seed);
    MtrNode * MakeZddTreeNode(unsigned int low, unsigned int size, unsigned int type);
    void zddPrintSubtable();
    void zddReduceHeap(Cudd_ReorderingType heuristic, int minsize);
    void zddShuffleHeap(int * permutation);
    void zddSymmProfile(int lower, int upper);
    void DumpDot(int n, ZDD* f, char ** inames, char ** onames, FILE * fp);

}; // Cudd


/**Class***********************************************************************

  Synopsis     [Class for BDD vectors.]

  Description  []

  SeeAlso      [BDD]

******************************************************************************/
class BDDvector {
    struct capsule {
	Cudd *manager;
	BDD *vect;
	int size;
	int ref;
    };
    capsule *p;
public:
    BDDvector(int size, Cudd *manager = 0, DdNode **nodes = 0);
    BDDvector(const BDDvector &from);
    ~BDDvector();
    BDDvector& operator=(const BDDvector& right);
    BDD& operator[](int i);
    int count() {return p->size;}
    Cudd *manager() {return p->manager;}
    void DumpDot(char ** inames = 0, char ** onames = 0, FILE * fp = stdout);
    void DumpDaVinci(
      char ** inames = 0,
      char ** onames = 0,
      FILE * fp = stdout);
    void DumpBlif(
      char ** inames = 0,
      char ** onames = 0,
      char * mname = 0,
      FILE * fp = stdout);
    void DumpDDcal(char ** inames = 0, char ** onames = 0, FILE * fp = stdout);
    void DumpFactoredForm(
      char ** inames = 0,
      char ** onames = 0,
      FILE * fp = stdout);
    BDD VectorSupport();
    int nodeCount();
    int VectorSupportSize();

}; // BDDvector


/**Class***********************************************************************

  Synopsis     [Class for ADD vectors.]

  Description  []

  SeeAlso      [ADD]

******************************************************************************/
class ADDvector {
    struct capsule {
	Cudd *manager;
	ADD *vect;
	int size;
	int ref;
    };
    capsule *p;
public:
    ADDvector(int size, Cudd *manager = 0, DdNode **nodes = 0);
    ADDvector(const ADDvector &from);
    ~ADDvector();
    ADDvector& operator=(const ADDvector& right);
    ADD& operator[](int i);
    int count() {return p->size;}
    Cudd *manager() {return p->manager;}
    void DumpDot(char ** inames = 0, char ** onames = 0, FILE * fp = stdout);
    void DumpDaVinci(
      char ** inames = 0,
      char ** onames = 0,
      FILE * fp = stdout);
    BDD VectorSupport();
    int VectorSupportSize();

}; // ADDvector


/**Class***********************************************************************

  Synopsis     [Class for ZDD vectors.]

  Description  []

  SeeAlso      [ZDD]

******************************************************************************/
class ZDDvector {
    struct capsule {
	Cudd *manager;
	ZDD *vect;
	int size;
	int ref;
    };
    capsule *p;
public:
    ZDDvector(int size, Cudd *manager = 0, DdNode **nodes = 0);
    ZDDvector(const ZDDvector &from);
    ~ZDDvector();
    ZDDvector& operator=(const ZDDvector& right);
    ZDD& operator[](int i);
    int count() {return p->size;}
    Cudd *manager() {return p->manager;}
    void DumpDot(char ** inames = 0, char ** onames = 0, FILE * fp = stdout);

}; // ZDDvector


extern void defaultError(string message);
extern int NextCube(DdGen * gen, int ** cube, CUDD_VALUE_TYPE * value);
extern int GenFree(DdGen * gen);
extern int IsGenEmpty(DdGen * gen);

#endif
