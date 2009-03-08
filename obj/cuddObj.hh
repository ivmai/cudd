/**CHeaderFile***************************************************************

  FileName    [cuddObj.hh]

  PackageName [cudd]

  Synopsis    [Class definitions for C++ object-oriented encapsulation of
  CUDD.]

  Description [Class definitions for C++ object-oriented encapsulation of
  CUDD.]

  SeeAlso     []

  Author      [Fabio Somenzi]

  Copyright   [Copyright (c) 1995-2004, Regents of the University of Colorado

  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions
  are met:

  Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

  Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

  Neither the name of the University of Colorado nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
  POSSIBILITY OF SUCH DAMAGE.]

  Revision    [$Id: cuddObj.hh,v 1.11 2009/02/21 19:41:38 fabio Exp fabio $]

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
    inline DdManager * checkSameManager(const DD &other) const;
    inline void checkReturnValue(const DdNode *result) const;
    inline void checkReturnValue(const int result, const int expected = 1)
	const;
public:
    DD(Cudd *ddManager, DdNode *ddNode);
    DD();
    DD(const DD &from);
    virtual ~DD();
    Cudd *manager() const;
    inline DdNode * getNode() const;
    int nodeCount() const;
    unsigned int NodeReadIndex() const;

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
    int operator==(const ABDD &other) const;
    int operator!=(const ABDD &other) const;
    void print(int nvars, int verbosity = 1) const;
    DdApaNumber ApaCountMinterm(int nvars, int * digits) const;
    void ApaPrintMinterm(int nvars, FILE * fp = stdout) const;
    void EpdPrintMinterm(int nvars, FILE * fp = stdout) const;
    BDD FindEssential() const;
    void PrintTwoLiteralClauses(char ** names, FILE * fp = stdout) const;
    BDD ShortestPath(int * weight, int * support, int * length) const;
    BDD LargestCube(int * length) const;
    int ShortestLength(int * weight) const;
    int EquivDC(const ABDD& G, const ABDD& D) const;
    double * CofMinterm() const;
    void PrintMinterm() const;
    double CountMinterm(int nvars) const;
    double CountPath() const;
    BDD Support() const;
    int SupportSize() const;
    void ClassifySupport(const ABDD& g, BDD* common, BDD* onlyF, BDD* onlyG)
	const;
    int CountLeaves() const;
    DdGen * FirstCube(int ** cube, CUDD_VALUE_TYPE * value) const;
    double Density(int nvars) const;

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
    int operator==(const BDD& other) const;
    int operator!=(const BDD& other) const;
    BDD operator=(const BDD& right);
    int operator<=(const BDD& other) const;
    int operator>=(const BDD& other) const;
    int operator<(const BDD& other) const;
    int operator>(const BDD& other) const;
    BDD operator!() const;
    BDD operator~() const;
    BDD operator*(const BDD& other) const;
    BDD operator*=(const BDD& other);
    BDD operator&(const BDD& other) const;
    BDD operator&=(const BDD& other);
    BDD operator+(const BDD& other) const;
    BDD operator+=(const BDD& other);
    BDD operator|(const BDD& other) const;
    BDD operator|=(const BDD& other);
    BDD operator^(const BDD& other) const;
    BDD operator^=(const BDD& other);
    BDD operator-(const BDD& other) const;
    BDD operator-=(const BDD& other);
    BDD AndAbstract(const BDD& g, const BDD& cube) const;
    BDD AndAbstractLimit(const BDD& g, const BDD& cube, unsigned int limit)
	const;
    BDD UnderApprox(
      int numVars,
      int threshold = 0,
      int safe = 0,
      double quality = 1.0) const;
    BDD OverApprox(
      int numVars,
      int threshold = 0,
      int safe = 0,
      double quality = 1.0) const;
    BDD RemapUnderApprox(int numVars, int threshold = 0, double quality = 1.0)
	const;
    BDD RemapOverApprox(int numVars, int threshold = 0, double quality = 1.0)
	const;
    BDD ExistAbstract(const BDD& cube) const;
    BDD XorExistAbstract(const BDD& g, const BDD& cube) const;
    BDD UnivAbstract(const BDD& cube) const;
    BDD BooleanDiff(int x) const;
    int VarIsDependent(const BDD& var) const;
    double Correlation(const BDD& g) const;
    double CorrelationWeights(const BDD& g, double * prob) const;
    BDD Ite(const BDD& g, const BDD& h) const;
    BDD IteConstant(const BDD& g, const BDD& h) const;
    BDD Intersect(const BDD& g) const;
    BDD And(const BDD& g) const;
    BDD AndLimit(const BDD& g, unsigned int limit) const;
    BDD Or(const BDD& g) const;
    BDD Nand(const BDD& g) const;
    BDD Nor(const BDD& g) const;
    BDD Xor(const BDD& g) const;
    BDD Xnor(const BDD& g) const;
    int Leq(const BDD& g) const;
    ADD Add() const;
    BDD Transfer(Cudd& destination) const;
    BDD ClippingAnd(const BDD& g, int maxDepth, int direction) const;
    BDD ClippingAndAbstract(const BDD& g, const BDD& cube, int maxDepth,
			    int direction) const;
    BDD Cofactor(const BDD& g) const;
    BDD Compose(const BDD& g, int v) const;
    BDD Permute(int * permut) const;
    BDD SwapVariables(BDDvector x, BDDvector y) const;
    BDD AdjPermuteX(BDDvector x) const;
    BDD VectorCompose(BDDvector vector) const;
    void ApproxConjDecomp(BDD* g, BDD* h) const;
    void ApproxDisjDecomp(BDD* g, BDD* h) const;
    void IterConjDecomp(BDD* g, BDD* h) const;
    void IterDisjDecomp(BDD* g, BDD* h) const;
    void GenConjDecomp(BDD* g, BDD* h) const;
    void GenDisjDecomp(BDD* g, BDD* h) const;
    void VarConjDecomp(BDD* g, BDD* h) const;
    void VarDisjDecomp(BDD* g, BDD* h) const;
    int IsVarEssential(int id, int phase) const;
    BDD Constrain(const BDD& c) const;
    BDD Restrict(const BDD& c) const;
    BDD NPAnd(const BDD& g) const;
    BDDvector ConstrainDecomp() const;
    BDDvector CharToVect() const;
    BDD LICompaction(const BDD& c) const;
    BDD Squeeze(const BDD& u) const;
    BDD Minimize(const BDD& c) const;
    BDD SubsetCompress(int nvars, int threshold) const;
    BDD SupersetCompress(int nvars, int threshold) const;
    BDD LiteralSetIntersection(const BDD& g) const;
    BDD PrioritySelect(BDDvector x, BDDvector y, BDDvector z, const BDD& Pi,
		       DD_PRFP Pifunc) const;
    BDD CProjection(const BDD& Y) const;
    int MinHammingDist(int *minterm, int upperBound) const;
    BDD Eval(int * inputs) const;
    BDD Decreasing(int i) const;
    BDD Increasing(int i) const;
    BDD SolveEqn(const BDD& Y, BDD* G, int ** yIndex, int n) const;
    BDD VerifySol(BDD* G, int * yIndex, int n) const;
    BDD SplitSet(BDDvector xVars, double m) const;
    BDD SubsetHeavyBranch(int numVars, int threshold) const;
    BDD SupersetHeavyBranch(int numVars, int threshold) const;
    BDD SubsetShortPaths(int numVars, int threshold, int hardlimit) const;
    BDD SupersetShortPaths(int numVars, int threshold, int hardlimit) const;
    void PrintCover() const;
    void PrintCover(const BDD& u) const;
    int EstimateCofactor(int i, int phase) const;
    int EstimateCofactorSimple(int i) const;
    void PickOneCube(char * string) const;
    BDD PickOneMinterm(BDDvector vars) const;
    DdGen * FirstNode(BDD* fnode) const;
    BDD zddIsop(const BDD& U, ZDD* zdd_I) const;
    BDD Isop(const BDD& U) const;
    ZDD PortToZdd() const;

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
    int operator==(const ADD& other) const;
    int operator!=(const ADD& other) const;
    ADD operator=(const ADD& right);
    // Relational operators
    int operator<=(const ADD& other) const;
    int operator>=(const ADD& other) const;
    int operator<(const ADD& other) const;
    int operator>(const ADD& other) const;
    // Arithmetic operators
    ADD operator-() const;
    ADD operator*(const ADD& other) const;
    ADD operator*=(const ADD& other);
    ADD operator+(const ADD& other) const;
    ADD operator+=(const ADD& other);
    ADD operator-(const ADD& other) const;
    ADD operator-=(const ADD& other);
    // Logical operators
    ADD operator~() const;
    ADD operator&(const ADD& other) const;
    ADD operator&=(const ADD& other);
    ADD operator|(const ADD& other) const;
    ADD operator|=(const ADD& other);
    ADD ExistAbstract(const ADD& cube) const;
    ADD UnivAbstract(const ADD& cube) const;
    ADD OrAbstract(const ADD& cube) const;
    ADD Plus(const ADD& g) const;
    ADD Times(const ADD& g) const;
    ADD Threshold(const ADD& g) const;
    ADD SetNZ(const ADD& g) const;
    ADD Divide(const ADD& g) const;
    ADD Minus(const ADD& g) const;
    ADD Minimum(const ADD& g) const;
    ADD Maximum(const ADD& g) const;
    ADD OneZeroMaximum(const ADD& g) const;
    ADD Diff(const ADD& g) const;
    ADD Agreement(const ADD& g) const;
    ADD Or(const ADD& g) const;
    ADD Nand(const ADD& g) const;
    ADD Nor(const ADD& g) const;
    ADD Xor(const ADD& g) const;
    ADD Xnor(const ADD& g) const;
    ADD Log() const;
    ADD FindMax() const;
    ADD FindMin() const;
    ADD IthBit(int bit) const;
    ADD ScalarInverse(const ADD& epsilon) const;
    ADD Ite(const ADD& g, const ADD& h) const;
    ADD IteConstant(const ADD& g, const ADD& h) const;
    ADD EvalConst(const ADD& g) const;
    int Leq(const ADD& g) const;
    ADD Cmpl() const;
    ADD Negate() const;
    ADD RoundOff(int N) const;
    BDD BddThreshold(CUDD_VALUE_TYPE value) const;
    BDD BddStrictThreshold(CUDD_VALUE_TYPE value) const;
    BDD BddInterval(CUDD_VALUE_TYPE lower, CUDD_VALUE_TYPE upper) const;
    BDD BddIthBit(int bit) const;
    BDD BddPattern() const;
    ADD Cofactor(const ADD& g) const;
    ADD Compose(const ADD& g, int v) const;
    ADD Permute(int * permut) const;
    ADD SwapVariables(ADDvector x, ADDvector y) const;
    ADD VectorCompose(ADDvector vector) const;
    ADD NonSimCompose(ADDvector vector) const;
    ADD Constrain(const ADD& c) const;
    ADD Restrict(const ADD& c) const;
    ADD MatrixMultiply(const ADD& B, ADDvector z) const;
    ADD TimesPlus(const ADD& B, ADDvector z) const;
    ADD Triangle(const ADD& g, ADDvector z) const;
    ADD Eval(int * inputs) const;
    int EqualSupNorm(const ADD& g, CUDD_VALUE_TYPE tolerance, int pr) const;

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
    int operator==(const ZDD& other) const;
    int operator!=(const ZDD& other) const;
    int operator<=(const ZDD& other) const;
    int operator>=(const ZDD& other) const;
    int operator<(const ZDD& other) const;
    int operator>(const ZDD& other) const;
    void print(int nvars, int verbosity = 1) const;
    ZDD operator*(const ZDD& other) const;
    ZDD operator*=(const ZDD& other);
    ZDD operator&(const ZDD& other) const;
    ZDD operator&=(const ZDD& other);
    ZDD operator+(const ZDD& other) const;
    ZDD operator+=(const ZDD& other);
    ZDD operator|(const ZDD& other) const;
    ZDD operator|=(const ZDD& other);
    ZDD operator-(const ZDD& other) const;
    ZDD operator-=(const ZDD& other);
    int Count() const;
    double CountDouble() const;
    ZDD Product(const ZDD& g) const;
    ZDD UnateProduct(const ZDD& g) const;
    ZDD WeakDiv(const ZDD& g) const;
    ZDD Divide(const ZDD& g) const;
    ZDD WeakDivF(const ZDD& g) const;
    ZDD DivideF(const ZDD& g) const;
    double CountMinterm(int path) const;
    BDD PortToBdd() const;
    ZDD Ite(const ZDD& g, const ZDD& h) const;
    ZDD Union(const ZDD& Q) const;
    ZDD Intersect(const ZDD& Q) const;
    ZDD Diff(const ZDD& Q) const;
    ZDD DiffConst(const ZDD& Q) const;
    ZDD Subset1(int var) const;
    ZDD Subset0(int var) const;
    ZDD Change(int var) const;
    void PrintMinterm() const;
    void PrintCover() const;

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
    PFC getHandler() const;
    DdManager *getManager() const {return p->manager;}
    inline void makeVerbose() {p->verbose = 1;}
    inline void makeTerse() {p->verbose = 0;}
    inline int isVerbose() const {return p->verbose;}
    inline void checkReturnValue(const DdNode *result) const;
    inline void checkReturnValue(const int result) const;
    Cudd& operator=(const Cudd& right);
    void info() const;
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
    int ReorderingStatus(Cudd_ReorderingType * method) const;
    void AutodynEnableZdd(Cudd_ReorderingType method);
    void AutodynDisableZdd();
    int ReorderingStatusZdd(Cudd_ReorderingType * method) const;
    int zddRealignmentEnabled() const;
    void zddRealignEnable();
    void zddRealignDisable();
    int bddRealignmentEnabled() const;
    void bddRealignEnable();
    void bddRealignDisable();
    ADD background();
    void SetBackground(ADD bg);
    unsigned int ReadCacheSlots() const;
    double ReadCacheUsedSlots() const;
    double ReadCacheLookUps() const;
    double ReadCacheHits() const;
    unsigned int ReadMinHit() const;
    void SetMinHit(unsigned int hr);
    unsigned int ReadLooseUpTo() const;
    void SetLooseUpTo(unsigned int lut);
    unsigned int ReadMaxCache() const;
    unsigned int ReadMaxCacheHard() const;
    void SetMaxCacheHard(unsigned int mc);
    int ReadSize() const;
    int ReadZddSize() const;
    unsigned int ReadSlots() const;
    unsigned int ReadKeys() const;
    unsigned int ReadDead() const;
    unsigned int ReadMinDead() const;
    int ReadReorderings() const;
    long ReadReorderingTime() const;
    int ReadGarbageCollections() const;
    long ReadGarbageCollectionTime() const;
    int ReadSiftMaxVar() const;
    void SetSiftMaxVar(int smv);
    int ReadSiftMaxSwap() const;
    void SetSiftMaxSwap(int sms);
    double ReadMaxGrowth() const;
    void SetMaxGrowth(double mg);
    MtrNode * ReadTree() const;
    void SetTree(MtrNode * tree);
    void FreeTree();
    MtrNode * ReadZddTree() const;
    void SetZddTree(MtrNode * tree);
    void FreeZddTree();
    int ReadPerm(int i) const;
    int ReadPermZdd(int i) const;
    int ReadInvPerm(int i) const;
    int ReadInvPermZdd(int i) const;
    BDD ReadVars(int i);
    CUDD_VALUE_TYPE ReadEpsilon() const;
    void SetEpsilon(CUDD_VALUE_TYPE ep);
    Cudd_AggregationType ReadGroupcheck() const;
    void SetGroupcheck(Cudd_AggregationType gc);
    int GarbageCollectionEnabled() const;
    void EnableGarbageCollection();
    void DisableGarbageCollection();
    int DeadAreCounted() const;
    void TurnOnCountDead();
    void TurnOffCountDead();
    int ReadRecomb() const;
    void SetRecomb(int recomb);
    int ReadSymmviolation() const;
    void SetSymmviolation(int symmviolation);
    int ReadArcviolation() const;
    void SetArcviolation(int arcviolation);
    int ReadPopulationSize() const;
    void SetPopulationSize(int populationSize);
    int ReadNumberXovers() const;
    void SetNumberXovers(int numberXovers);
    unsigned long ReadMemoryInUse() const;
    long ReadPeakNodeCount() const;
    long ReadNodeCount() const;
    long zddReadNodeCount() const;
    void AddHook(DD_HFP f, Cudd_HookType where);
    void RemoveHook(DD_HFP f, Cudd_HookType where);
    int IsInHook(DD_HFP f, Cudd_HookType where) const;
    void EnableReorderingReporting();
    void DisableReorderingReporting();
    int ReorderingReporting();
    int ReadErrorCode() const;
    void ClearErrorCode();
    FILE *ReadStdout() const;
    void SetStdout(FILE *);
    FILE *ReadStderr() const;
    void SetStderr(FILE *);
    unsigned int ReadNextReordering() const;
    double ReadSwapSteps() const;
    unsigned int ReadMaxLive() const;
    void SetMaxLive(unsigned int);
    unsigned long ReadMaxMemory() const;
    void SetMaxMemory(unsigned long);
    int bddBindVar(int);
    int bddUnbindVar(int);
    int bddVarIsBound(int) const;
    ADD Walsh(ADDvector x, ADDvector y);
    ADD addResidue(int n, int m, int options, int top);
    int ApaNumberOfDigits(int binaryDigits) const;
    DdApaNumber NewApaNumber(int digits) const;
    void ApaCopy(int digits, DdApaNumber source, DdApaNumber dest) const;
    DdApaDigit ApaAdd(int digits, DdApaNumber a, DdApaNumber b, DdApaNumber
		      sum) const;
    DdApaDigit ApaSubtract(int digits, DdApaNumber a, DdApaNumber b,
			   DdApaNumber diff) const;
    DdApaDigit ApaShortDivision(int digits, DdApaNumber dividend, DdApaDigit
				divisor, DdApaNumber quotient) const;
    void ApaShiftRight(int digits, DdApaDigit in, DdApaNumber a, DdApaNumber
		       b) const;
    void ApaSetToLiteral(int digits, DdApaNumber number, DdApaDigit literal)
      const;
    void ApaPowerOfTwo(int digits, DdApaNumber number, int power) const;
    void ApaPrintHex(FILE * fp, int digits, DdApaNumber number) const;
    void ApaPrintDecimal(FILE * fp, int digits, DdApaNumber number) const;
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
    BDD Inequality(int c, BDDvector x, BDDvector y);
    BDD Disequality(int c, BDDvector x, BDDvector y);
    BDD Interval(BDDvector x, unsigned int lowerB, unsigned int upperB);
    ADD Hamming(ADDvector xVars, ADDvector yVars);
    // void Read(FILE * fp, ADD* E, ADD** x, ADD** y, ADD** xn, ADD** yn_, int * nx, int * ny, int * m, int * n, int bx, int sx, int by, int sy);
    // void Read(FILE * fp, BDD* E, BDD** x, BDD** y, int * nx, int * ny, int * m, int * n, int bx, int sx, int by, int sy);
    void ReduceHeap(Cudd_ReorderingType heuristic, int minsize);
    void ShuffleHeap(int * permutation);
    void SymmProfile(int lower, int upper) const;
    unsigned int Prime(unsigned int pr) const;
    int SharingSize(DD* nodes, int n) const;
    BDD bddComputeCube(BDD * vars, int * phase, int n);
    ADD addComputeCube(ADD * vars, int * phase, int n);
    int NextNode(DdGen * gen, BDD * nnode);
    BDD IndicesToCube(int * array, int n);
    void PrintVersion(FILE * fp) const;
    double AverageDistance() const;
    long Random();
    void Srandom(long seed);
    MtrNode * MakeZddTreeNode(unsigned int low, unsigned int size, unsigned int type);
    void zddPrintSubtable() const;
    void zddReduceHeap(Cudd_ReorderingType heuristic, int minsize);
    void zddShuffleHeap(int * permutation);
    void zddSymmProfile(int lower, int upper) const;
  //void DumpDot(int n, ZDD* f, char ** inames, char ** onames, FILE * fp);

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
    BDD& operator[](int i) const;
    int count() const {return p->size;}
    Cudd *manager() const {return p->manager;}
    void DumpDot(char ** inames = 0, char ** onames = 0, FILE * fp = stdout)
	const;
    void DumpDaVinci(
      char ** inames = 0,
      char ** onames = 0,
      FILE * fp = stdout) const;
    void DumpBlif(
      char ** inames = 0,
      char ** onames = 0,
      char * mname = 0,
      FILE * fp = stdout,
      int mv = 0) const;
    void DumpDDcal(char ** inames = 0, char ** onames = 0, FILE * fp = stdout)
	const;
    void DumpFactoredForm(
      char ** inames = 0,
      char ** onames = 0,
      FILE * fp = stdout) const;
    BDD VectorSupport() const;
    int nodeCount() const;
    int VectorSupportSize() const;

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
    ADD& operator[](int i) const;
    int count() const {return p->size;}
    Cudd *manager() const {return p->manager;}
    void DumpDot(char ** inames = 0, char ** onames = 0, FILE * fp = stdout)
	const;
    void DumpDaVinci(
      char ** inames = 0,
      char ** onames = 0,
      FILE * fp = stdout) const;
    BDD VectorSupport() const;
    int VectorSupportSize() const;

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
    ZDD& operator[](int i) const;
    int count() const {return p->size;}
    Cudd *manager() const {return p->manager;}
    void DumpDot(char ** inames = 0, char ** onames = 0, FILE * fp = stdout)
	const;

}; // ZDDvector


extern void defaultError(string message);
extern int NextCube(DdGen * gen, int ** cube, CUDD_VALUE_TYPE * value);
extern int GenFree(DdGen * gen);
extern int IsGenEmpty(DdGen * gen);

#endif
