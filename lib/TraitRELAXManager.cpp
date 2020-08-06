#include "TraitRELAXManager.h"

// constructor
TraitRELAXManager::TraitRELAXManager(BppApplication* parameters) :
  traitRELAXParameters_(),
  traitRELAXLikelihoodFunction_(),
  nullLogl_(0),
  alternativeLogl_(0),
  gCode_(),
  binaryAlphabet_(),
  characterData_(),
  tree_(),
  codonAlphabet_(),
  sequenceData_(),
  characterModel_(),
  sequenceModel_(),
  rDist_()
{
  traitRELAXParameters_ = parameters;
}

// destructor
TraitRELAXManager::~TraitRELAXManager()
{
  traitRELAXParameters_->done();
  if (gCode_) delete gCode_;          
  if (binaryAlphabet_) delete binaryAlphabet_;
  if (characterData_) delete characterData_;
  if (tree_) delete tree_;
  if (codonAlphabet_->getNucleicAlphabet()) delete codonAlphabet_->getNucleicAlphabet();
  if (codonAlphabet_) delete codonAlphabet_;
  if (sequenceData_) delete sequenceData_;
  if (characterModel_) delete characterModel_;
  if (sequenceModel_) delete sequenceModel_;
  if (rDist_) delete rDist_;
  if (traitRELAXLikelihoodFunction_) delete traitRELAXLikelihoodFunction_;
}

/******************************************************************************/
/**************************** TraitRELAX functions ****************************/
/******************************************************************************/

// creates a codon alphabet
void TraitRELAXManager::setCodonAlphabet()
{
  map<string, string> alphabetSettings;
  alphabetSettings["alphabet"] = "Codon(letter=DNA)";
  const Alphabet* alphabet = SequenceApplicationTools::getAlphabet(alphabetSettings, "", false);
  codonAlphabet_ = dynamic_cast<const CodonAlphabet *>(alphabet);
}

/******************************************************************************/

// processes input character trait data
void TraitRELAXManager::setCharacterData()
{
  binaryAlphabet_ = new BinaryAlphabet();
  string charDataFilePath = ApplicationTools::getAFilePath("input.character.file", traitRELAXParameters_->getParams(), true, true, "", true, "none", 1);
  string sequenceFormat = ApplicationTools::getStringParameter("input.character.format", traitRELAXParameters_->getParams(), "Fasta()", "", true, 1);
  BppOAlignmentReaderFormat bppoReader(1);
  unique_ptr<IAlignment> iAln(bppoReader.read(sequenceFormat));
  map<string, string> args(bppoReader.getUnparsedArguments());
  ApplicationTools::displayResult("character data file ", charDataFilePath);
  ApplicationTools::displayResult("chatacter data format ", iAln->getFormatName());
  const SequenceContainer *charCont = iAln->readAlignment(charDataFilePath, binaryAlphabet_);
  VectorSiteContainer *sites = new VectorSiteContainer(*dynamic_cast<const OrderedSequenceContainer *>(charCont));
  delete charCont;
  characterData_ = sites;
}

/******************************************************************************/
/************************** Environment setup *********************************/
/******************************************************************************/

// processes input sequence data
void TraitRELAXManager::setSequenceData()
{
  VectorSiteContainer *allSites = SequenceApplicationTools::getSiteContainer(codonAlphabet_, traitRELAXParameters_->getParams());            // here, gaps will be converted to unknown characters
  VectorSiteContainer *sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, traitRELAXParameters_->getParams(), "", true, false); // convert the alignemnt to condensed format of unique sites
  delete allSites;                                                                                                          // delete the non-condenced intance of the sequence data
  SiteContainerTools::changeGapsToUnknownCharacters(*sites);                                                                // convert gaps to unknown characters (as done in bppML.cpp)
  sequenceData_ = sites;
}

/******************************************************************************/

// gives names to internal nodes for the trait histories internal nodes identification
void TraitRELAXManager::giveNamesToInternalNodes(Tree *tree)
{
  TreeTemplate<Node> *ttree = dynamic_cast<TreeTemplate<Node> *>(tree);
  vector<Node *> nodes = ttree->getNodes();
  for (size_t i = 0; i < nodes.size(); ++i)
  {
    if (!nodes[i]->hasName())
      nodes[i]->setName("_baseInternal_" + TextTools::toString(nodes[i]->getId()));
  }
}

/******************************************************************************/

// processes input tree
void TraitRELAXManager::setTree()
{
  Tree *tree = PhylogeneticsApplicationTools::getTree(traitRELAXParameters_->getParams());
  giveNamesToInternalNodes(tree);
  tree_ = tree;
}

/******************************************************************************/

// creates the character trait model
void TraitRELAXManager::setCharacterModel(DRTreeParsimonyScore *mpData)
{
  // create the model
  SubstitutionModel *model = new TwoParameterBinarySubstitutionModel(binaryAlphabet_);

  // compute the maximum parsimony score and set the lower and upper bounds on mu (=rate) as mp/tree_size, 2*mp/tree_size
  VDouble treeBranches = dynamic_cast<TreeTemplate<Node> *>(tree_)->getBranchLengths();
  double treeSize = 0;
  for (size_t i = 0; i < treeBranches.size(); ++i)
  {
    treeSize += treeBranches[i];
  }
  double characterMuLb = mpData->getScore() / treeSize;
  double characterMuUb = 4 * characterMuLb; // used 4*lb instead of 2*lb because factor of 2 is not enough (optimization converges to upper bound)

  // set the initial values of the model
  double mu = (characterMuLb + characterMuUb) / 2;
  if (!ApplicationTools::getBooleanParameter("character_model.set_initial_parameters", traitRELAXParameters_->getParams(), true, "", true, false))
  {
    // set the value of mu to be the middle of the interval
    model->setParameterValue(string("mu"), mu);
	// estimate the initial frequencies as observedPseudoCount with pseudocount as 1 to avoid possible case of frequency = 0
    model->setFreqFromData(dynamic_cast<const SequenceContainer &>(*characterData_), 1); // the second arguemnt stands for pesudocount 1
  }
  else
  {
	double init_mu = ApplicationTools::getDoubleParameter("character_model.mu", traitRELAXParameters_->getParams(), 1);
    double init_pi0 = ApplicationTools::getDoubleParameter("character_model.pi0", traitRELAXParameters_->getParams(), 0.5);
    model->setParameterValue(string("mu"), init_mu);
	model->setParameterValue(string("pi0"), init_pi0);
	mu = init_mu;
  }
  if (mu < characterMuLb)
  {
    dynamic_cast<TwoParameterBinarySubstitutionModel *>(model)->setMuBounds(mu - 0.001, characterMuUb);
  }
  else if (mu > characterMuUb)
  {
    dynamic_cast<TwoParameterBinarySubstitutionModel *>(model)->setMuBounds(characterMuLb, mu + 0.001);
  }
  else
  {
	dynamic_cast<TwoParameterBinarySubstitutionModel *>(model)->setMuBounds(characterMuLb, characterMuUb);
  }
  
  characterModel_ = dynamic_cast<TransitionModel *>(model);
}

/******************************************************************************/

// sets initial partition of the branches as input of the sequence model based on a maximum parsimony solution
void TraitRELAXManager::setMpPartition(DRTreeParsimonyScore *mpData)
{
  mpData->computeSolution();
  const Tree &solution = mpData->getTree();
  vector<const Node *> nodes = (dynamic_cast<const TreeTemplate<Node> &>(solution)).getNodes();

  // set assignment to branches
  string character0NodesIds, character1NodesIds = "";
  for (size_t i = 0; i < nodes.size(); ++i)
  {
    if (!tree_->isRoot(nodes[i]->getId()))
    {
      int nodeState = dynamic_cast<const BppInteger *>(nodes[i]->getNodeProperty("state"))->getValue();
      if (nodeState == 0)
      {
        character0NodesIds = character0NodesIds + TextTools::toString(nodes[i]->getId()) + ",";
      }
      else
      {
        character1NodesIds = character1NodesIds + TextTools::toString(nodes[i]->getId()) + ",";
      }
    }
  }
  traitRELAXParameters_->getParam("model1.nodes_id") = character0NodesIds.substr(0, character0NodesIds.length() - 1);
  traitRELAXParameters_->getParam("model2.nodes_id") = character1NodesIds.substr(0, character1NodesIds.length() - 1);
}

/******************************************************************************/

// creates the sequence model
void TraitRELAXManager::setSequenceModel(DRTreeParsimonyScore *mpData)
{
  // set initial partition, based on maximum parsimony
  setMpPartition(mpData); // the partition is set on tree

  MixedSubstitutionModelSet *modelSet; 
  // TransitionModel *model; // commented out until further consult with Itay as to weather this should be implemented

  if (!ApplicationTools::getBooleanParameter("sequence_model.set_initial_parameters", traitRELAXParameters_->getParams(), true, "", true, false))
  {
    traitRELAXParameters_->getParam("model1") = "RELAX(kappa=1,p=0.1,omega1=1.0,omega2=2.0,theta1=0.5,theta2=0.8,k=1,frequencies=F3X4,initFreqs=observed,initFreqs.observedPseudoCount=1)";
    traitRELAXParameters_->getParam("model2") = "RELAX(kappa=RELAX.kappa_1,p=RELAX.p_1,omega1=RELAX.omega1_1,omega2=RELAX.omega2_1,theta1=RELAX.theta1_1,theta2=RELAX.theta2_1,k=1,1_Full.theta=RELAX.1_Full.theta_1,1_Full.theta1=RELAX.1_Full.theta1_1,1_Full.theta2=RELAX.1_Full.theta2_1,2_Full.theta=RELAX.2_Full.theta_1,2_Full.theta1=RELAX.2_Full.theta1_1,2_Full.theta2=RELAX.2_Full.theta2_1,3_Full.theta=RELAX.3_Full.theta_1,3_Full.theta1=RELAX.3_Full.theta1_1,3_Full.theta2=RELAX.3_Full.theta2_1,frequencies=F3X4,initFreqs=observed,initFreqs.observedPseudoCount=1)";
  }
  traitRELAXParameters_->getParam("nonhomogeneous") = "general";
  traitRELAXParameters_->getParam("nonhomogeneous.number_of_models") = "2";
  traitRELAXParameters_->getParam("nonhomogeneous.stationarity") = "yes"; 					 // constrain root frequencies to be the same as stationary (since RELAX is a time reversible model, this should not cause issues)


	int LLApproach = ApplicationTools::getIntParameter("likelihood.computation.approach", traitRELAXParameters_->getParams(), 3);

  // commented out until further consult with Itay as to weather this should be implemented
  // likelihood approach 1 (RELAX HyPhy default)- set likelihood computation to enable shifts in the selective regime in different branches in the phylogeny with respect to each site
  if (LLApproach == 1)
  {
    throw Exception("Random effects likelihood approach is currently disabled. PLease choose 2 or 3.");
  }
  // likelihood approach 3 (TraitRELAX default) - set likelihood computation to restrict the same selective regime for each branch along the phylogeny with respect to each site
  if (LLApproach == 3)
  {
    traitRELAXParameters_->getParam("site.number_of_paths") = "2";                               // the 3rd path mapping omega3 in the branches under chatacter states 0 and 1 is imlied by the other two paths
    traitRELAXParameters_->getParam("site.path1") = "model1[YN98.omega_1]&model2[YN98.omega_1]"; // map omega1 in the branches under character state 0 (=model1) to omega1 in the branches under character state 1 (=model2)
    traitRELAXParameters_->getParam("site.path2") = "model1[YN98.omega_2]&model2[YN98.omega_2]"; // do the same for omega2
  }
  // likelihood approach 2 (Bio++ default, no need to set) - set likelihood computation to enable shifts in the selective regime between branches in the phylogeny with respect to each site only upon transition from FG to BG category

  // create the set of models
  string codeDesc = ApplicationTools::getStringParameter("genetic_code", traitRELAXParameters_->getParams(), "Standard", "", true, true);
  gCode_ = SequenceApplicationTools::getGeneticCode(codonAlphabet_->getNucleicAlphabet(), codeDesc);
  modelSet = dynamic_cast<MixedSubstitutionModelSet *>(PhylogeneticsApplicationTools::getSubstitutionModelSet(codonAlphabet_, gCode_, sequenceData_, traitRELAXParameters_->getParams()));
  sequenceModel_ = modelSet;
}

/******************************************************************************/

// setting up environment for traitrelax execution
void TraitRELAXManager::init()
{
    // process seed from parameter file, if exists
    double seed;
    seed = ApplicationTools::getDoubleParameter("seed", traitRELAXParameters_->getParams(), 1);
    if (seed == 1)
    {
      // else, choose a random seed
      seed = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0);
    }
    cout << "seed=" << seed << endl;
    RandomTools::setSeed(static_cast<long>(seed));

    /* process character data */
    setCharacterData();

    /* process tree */
    setTree();

    /* process codon alignment */
    setCodonAlphabet();
    setSequenceData();

    /* compute the maximum parsimony  for the purpose of setting bounds on the rate parameter of the character model and an intial tree partition for the starting point */
    DRTreeParsimonyScore* mpData = new DRTreeParsimonyScore(*tree_, dynamic_cast<const SiteContainer &>(*characterData_));

    /* set the character model */
    setCharacterModel(mpData);

    /* set the sequence model */
    setSequenceModel(mpData);
    delete mpData;

    /* set the joint likelihood function instance */
    rDist_ = new ConstantRateDistribution();

    bool debug = ApplicationTools::getBooleanParameter("joint.likelihood.debug", traitRELAXParameters_->getParams(), false, "", true, 1);
    traitRELAXLikelihoodFunction_ = new JointLikelihoodFunction(traitRELAXParameters_, tree_, characterData_, characterModel_, sequenceData_, sequenceModel_, rDist_, debug);
}

/******************************************************************************/

void TraitRELAXManager::reportParameters(map <string, double> parameterValues)
{
    // report the optimal log likelihood and parameters of the alternative model
  for (map<string, double>::iterator it = parameterValues.begin(); it != parameterValues.end(); it++)
  {
    if ((it->first.find("Log likelihood") == std::string::npos))
    {
      ApplicationTools::displayResult(it->first, TextTools::toString(it->second, 15));
    }
  }
	// now report the log likelihood value
	cout << "\n" << endl;
	ApplicationTools::displayResult("Character Log likelihood", TextTools::toString(-parameterValues["Character Log likelihood"], 15));
	ApplicationTools::displayResult("Sequence Log likelihood", TextTools::toString(-parameterValues["Sequence Log likelihood"], 15));
	ApplicationTools::displayResult("Overall Log likelihood", TextTools::toString(-parameterValues["Overall Log likelihood"], 15));
}

/******************************************************************************/
/************************** Null model optimization ***************************/
/******************************************************************************/

map<string, double> TraitRELAXManager::optimizeNullModel()
{ 
    cout << "\n*********************************** Null model fitting **********************************" << endl;
    unsigned int verbose = static_cast<unsigned int>(ApplicationTools::getIntParameter("optimization.verbose", traitRELAXParameters_->getParams(), 0));
    bool report = false;
    if (verbose)
    {
      cout << "\n**** Initial parameters ****" << endl;
      report = true;
      traitRELAXLikelihoodFunction_->getModelParameters(report);
    }
    ParameterList emptyParametersList;

    /* fit the null model: separate optimization of the character model and the sequence model, where the selection intensity parameter k is 1 */
    traitRELAXParameters_->startTimer();
    traitRELAXLikelihoodFunction_->setHypothesis(JointLikelihoodFunction::Hypothesis(0));
    traitRELAXLikelihoodFunction_->setOptimizationScope(JointLikelihoodFunction::OptimizationScope(3));
    traitRELAXLikelihoodFunction_->fireParameterChanged(emptyParametersList);
    nullLogl_ = -traitRELAXLikelihoodFunction_->getValue();
    map<string, double> bestModelParameters = traitRELAXLikelihoodFunction_->getModelParameters(false);
    reportParameters(bestModelParameters);
    traitRELAXParameters_->done();
    return bestModelParameters;
}

/******************************************************************************/
/********************** Alternative model optimization ************************/
/******************************************************************************/

// initalizes grid for character model optimization under the alternative hypothesis
void TraitRELAXManager::initCharacterGrid(VVDouble &grid, double pi0Lb, double pi0Ub, double muLb, double muUb, uint gridSize)
{
  double stepSize;
  double next;

  /* define a set of possible assignments for rate */
  if (muLb < 0 && muUb > 1) // in case of a range that allows both relative character substitution rate lower and greater than 1, set the jumps in the grid to be uneven in order to sample same amount of rates < 1 and rates > 1
  {
    double logUb = log(muUb);
    double logLb = log(muLb);
    stepSize = (logUb - logLb) / (gridSize - 1);
    next = logLb;
    for (size_t i = 0; i < gridSize; ++i)
    {
      grid[0][i] = exp(next);
      next += stepSize;
    }
  }
  else
  {
    stepSize = (muUb - muLb) / (gridSize - 1);
    next = muLb;
    for (size_t j = 0; j < gridSize; ++j)
    {
      grid[0][j] = next;
      next += stepSize;
    }
  }

  /* define a set of possible assignments for kappa */
  stepSize = (pi0Ub - pi0Lb) / (gridSize - 1);
  next = pi0Lb;
  for (size_t k = 0; k < gridSize; ++k)
  {
    grid[1][k] = next;
    next += stepSize;
  }
}

/******************************************************************************/

// obtains neighbors of a grid visited point
VVDouble TraitRELAXManager::getNeighbors(VVDouble &grid, double currentMu, double currentPi0)
{
  VVDouble neighbors;
  neighbors.clear();
  neighbors.resize(2, VDouble(3));
  // get the location of the current point in the grid
  uint muPos = 0;
  uint pi0Pos = 0;
  for (uint i = 0; i < grid.size(); i++)
  {
    if (grid[0][i] == currentMu)
    {
      muPos = i;
      break;
    }
  }
  for (uint j = 0; j < grid.size(); j++)
  {
    if (grid[1][j] == currentPi0)
    {
      pi0Pos = j;
      break;
    }
  }
  // return a double-list of all its surronding points
  int muMin, muMax, pi0Min, pi0Max;
  muMin = muPos - 1;
  muMax = muPos + 1;
  if (muPos == 0)
  {
    muMin = muMin;
  }
  if (muPos == grid.size())
  {
    muMax = muPos;
  }
  pi0Min = pi0Pos - 1;
  pi0Max = pi0Pos + 1;
  if (pi0Pos == 0)
  {
    pi0Min = pi0Min;
  }
  if (pi0Pos == grid.size())
  {
    pi0Max = pi0Pos;
  }
  for (uint i = muMin; i < static_cast<uint>(muMax); i++)
  {
    for (uint j = pi0Min; j < static_cast<uint>(pi0Max); j++)
    {
      neighbors[0][i - muMin] = grid[0][i];
      neighbors[1][j - pi0Min] = grid[1][j];
    }
  }
  return neighbors;
}

/******************************************************************************/

// optimizes the character model under the alternative hypothesis (in which the sequence evolution is linked to the character evolution of the trait) using one dimensional grid optimization
void TraitRELAXManager::optimizeAlternativeCharacterModelByGrid(map<string, double> &optimalValues, uint verbose, uint gridSize, bool firstCycle)
{
  cout << "* Optimizng joint likelihood function with respect to character parameters using grid *\n" << endl;

  // set the grid bounds on the rate and kappa parameters of the character model
  //const IntervalConstraint *muBounds = dynamic_cast<const IntervalConstraint *>(characterModel->getParameter("mu").getConstraint());
  std::shared_ptr<IntervalConstraint> muBounds = std::dynamic_pointer_cast<IntervalConstraint>(characterModel_->getParameter("mu").getConstraint());
  double muLb = muBounds->getLowerBound() + 0.0001;
  double muUb = muBounds->getUpperBound() - 0.0001;
  double pi0Lb = 0.1;
  double pi0Ub = 0.9; // altered from 999 to avoid radical values (-> failure to simulate history along a branch) set the upper bound on kappa such that pi0 = 1/(kappa+1) in [0.001,0.999] -> kappa in [0,999] (exclude 0 and 1 from the bouds to avoid absorving states in the character model)
  VVDouble grid;
  grid.clear();
  grid.resize(2, VDouble(gridSize));
  initCharacterGrid(grid, pi0Lb, pi0Ub, muLb, muUb, gridSize);
  
  // scan the grid and for each assignment of (rate, kappa), compute the joint likelihood and update the pair yieldling the best likelihood accordingly
  Parameter mu = traitRELAXLikelihoodFunction_->getParameter("TwoParameterBinary.mu");
  Parameter pi0 = traitRELAXLikelihoodFunction_->getParameter("TwoParameterBinary.pi0");
  double currentMu = mu.getValue();
  double currentPi0 = pi0.getValue();
  double currentLogL = -traitRELAXLikelihoodFunction_->getValue();
  optimalValues["mu"] = currentMu;
  optimalValues["pi0"] = currentPi0;
  optimalValues["log_likelihood"] = currentLogL;
  if (firstCycle)
  {
    for (size_t i = 0; i < gridSize; ++i)
    {
      currentMu = grid[0][i];
      mu.setValue(currentMu); // set the value of the rate
      for (size_t j = 0; j < gridSize; ++j)
      {
        currentPi0 = grid[1][j];
        pi0.setValue(currentPi0); // set the value of kappa
        ParameterList paramsToUpdate;
        paramsToUpdate.addParameter(mu);
        paramsToUpdate.addParameter(pi0);
        traitRELAXLikelihoodFunction_->setParametersValues(paramsToUpdate); // wrong values of parameters are set here. trigger likelihood computation
        currentLogL = -traitRELAXLikelihoodFunction_->getValue();
        if (currentLogL > optimalValues["log_likelihood"])
        {
          optimalValues["log_likelihood"] = currentLogL;
          optimalValues["mu"] = currentMu;
          optimalValues["pi0"] = currentPi0;
        }
      }
    }
  }
  else // i case of a non initial cycle, only traverse the neighbors of the current point in the likelihood surface
  {
    VVDouble neighbors = getNeighbors(grid, currentMu, currentPi0);
    for (uint i = 0; i < neighbors.size(); i++)
    {
      currentMu = grid[0][i];
      mu.setValue(currentMu); // set the value of the rate
      currentPi0 = grid[1][i];
      pi0.setValue(currentPi0); // set the value of kappa
      ParameterList paramsToUpdate;
      paramsToUpdate.addParameter(mu);
      paramsToUpdate.addParameter(pi0);
      traitRELAXLikelihoodFunction_->setParametersValues(paramsToUpdate); // wrong values of parameters are set here. trigger likelihood computation
      currentLogL = -traitRELAXLikelihoodFunction_->getValue();
      if (currentLogL > optimalValues["log_likelihood"])
      {
        optimalValues["log_likelihood"] = currentLogL;
        optimalValues["mu"] = currentMu;
        optimalValues["pi0"] = currentPi0;
      }
    }
  }
  // set the optimial values to the model
  traitRELAXLikelihoodFunction_->setParameterValue("TwoParameterBinary.mu", optimalValues["mu"]);
  traitRELAXLikelihoodFunction_->setParameterValue("TwoParameterBinary.pi0", optimalValues["pi0"]);
}

/******************************************************************************/

// optimizes the character model under the alternative hypothesis (in which the sequence evolution is linked to the character evolution of the trait) using one dimensional brent optimization
void TraitRELAXManager::optimizeAlternativeCharacterModelByOneDimBrent(map<string, double> &optimalValues, uint verbose)
{
  cout << "* Optimizing joint likelihood function with respect to character parameters using one dimentional brent *\n" << endl;

  // set the brent one dimontional optimizer
  BrentOneDimension *characterParametersOptimizer = new BrentOneDimension(traitRELAXLikelihoodFunction_);
  characterParametersOptimizer->setBracketing(BrentOneDimension::BRACKET_INWARD);
  characterParametersOptimizer->getStopCondition()->setTolerance(0.01);               // set the tolerance to be slighly less strict to account for the instability of the joint likelihood function
  characterParametersOptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO); // 3.4.19 - keren - OTHERWISE, CONSTRAINT IS EXCEEDED
  characterParametersOptimizer->setProfiler(0);
  characterParametersOptimizer->setMessageHandler(0);
  characterParametersOptimizer->setVerbose(1);
  ParameterList pi0, mu;
  pi0.addParameter(traitRELAXLikelihoodFunction_->getParameter("TwoParameterBinary.pi0"));
  //const IntervalConstraint *pi0Bounds = dynamic_cast<const IntervalConstraint *>(traitRELAXLikelihoodFunction_->getParameter("TwoParameterBinary.pi0").getConstraint());
  std::shared_ptr<IntervalConstraint> pi0Bounds = std::dynamic_pointer_cast<IntervalConstraint>(traitRELAXLikelihoodFunction_->getParameter("TwoParameterBinary.pi0").getConstraint());
  mu.addParameter(traitRELAXLikelihoodFunction_->getParameter("TwoParameterBinary.mu"));
  //const IntervalConstraint *muBounds = dynamic_cast<const IntervalConstraint *>(traitRELAXLikelihoodFunction_->getParameter("TwoParameterBinary.mu").getConstraint());
  std::shared_ptr<IntervalConstraint> muBounds = std::dynamic_pointer_cast<IntervalConstraint>(traitRELAXLikelihoodFunction_->getParameter("TwoParameterBinary.mu").getConstraint());
  double prevLogLikelihood, currLogLikelihood;
  do
  {
	  prevLogLikelihood = -traitRELAXLikelihoodFunction_->getValue();
	  
	  // optimize the joint model with respect to pi0
	  characterParametersOptimizer->setInitialInterval(pi0Bounds->getLowerBound(), pi0Bounds->getUpperBound()); // search within stricter bounds that the actual ones of pi0 to avoid failute of stochasitc mapping
	  characterParametersOptimizer->init(pi0);
	  characterParametersOptimizer->optimize();

	  optimalValues["pi0"] = traitRELAXLikelihoodFunction_->getParameter("TwoParameterBinary.pi0").getValue();

	  // optimize the model with respect to the mu
	  characterParametersOptimizer->setInitialInterval(muBounds->getLowerBound() + 0.0001, muBounds->getUpperBound() - 0.0001);
	  characterParametersOptimizer->init(mu);
	  characterParametersOptimizer->optimize();
	  optimalValues["mu"] = traitRELAXLikelihoodFunction_->getParameter("TwoParameterBinary.mu").getValue();
	  optimalValues["log_likelihood"] = traitRELAXLikelihoodFunction_->getCharacterLikelihoodFunction()->getValue();
	  
	  currLogLikelihood = -traitRELAXLikelihoodFunction_->getValue();
  } while (abs(currLogLikelihood - prevLogLikelihood) > 0.000001);
  delete characterParametersOptimizer;
}

/******************************************************************************/

// optimizes the character model under the alternative hypothesis (in which the sequence evolution is linked to the character evolution of the trait) using one dimensional powell optimization
void TraitRELAXManager::optimizeAlternativeCharacterModelByPowell(map<string, double> &optimalValues, uint verbose)
{
  cout << "* Optimizing joint likelihood function with respect to character parameters using two dimentional brent (i.e., powell) *\n" << endl;
  // set the brent two dimontional optimizer
  PowellMultiDimensions *characterParametersOptimizer = new PowellMultiDimensions(traitRELAXLikelihoodFunction_);
  ParameterList parametersToEstimate;
  parametersToEstimate.addParameter(traitRELAXLikelihoodFunction_->getParameter("TwoParameterBinary.mu"));
  parametersToEstimate.addParameter(traitRELAXLikelihoodFunction_->getParameter("TwoParameterBinary.pi0"));
  characterParametersOptimizer->setProfiler(0);
  characterParametersOptimizer->setMessageHandler(0);
  characterParametersOptimizer->setMaximumNumberOfEvaluations(1000);
  characterParametersOptimizer->getStopCondition()->setTolerance(0.01);
  characterParametersOptimizer->setVerbose(1);
  characterParametersOptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  characterParametersOptimizer->init(parametersToEstimate);

  double prevLogLikelihood = -traitRELAXLikelihoodFunction_->getValue();
  double currLogLikelihood = -traitRELAXLikelihoodFunction_->getValue();
  size_t index = 1;
  do
  {
    cout << "Optimization cycle: " << TextTools::toString(index) << endl;
    index = index + 1;
    prevLogLikelihood = -traitRELAXLikelihoodFunction_->getValue();


    characterParametersOptimizer->optimize();

    optimalValues["pi0"] = traitRELAXLikelihoodFunction_->getParameter("TwoParameterBinary.pi0").getValue();
    optimalValues["mu"] = traitRELAXLikelihoodFunction_->getParameter("TwoParameterBinary.mu").getValue();
    optimalValues["log_likelihood"] = traitRELAXLikelihoodFunction_->getCharacterLikelihoodFunction()->getValue();

    currLogLikelihood = -traitRELAXLikelihoodFunction_->getValue();
    ApplicationTools::displayResult("Current log likelihood", TextTools::toString(-currLogLikelihood, 15));
    ApplicationTools::displayResult("Current diff", TextTools::toString((currLogLikelihood-prevLogLikelihood), 15));

  } while (currLogLikelihood - prevLogLikelihood > 0.01);
  delete characterParametersOptimizer;

}

/******************************************************************************/

// optimizes the character model under the alternative hypothesis (in which the sequence evolution is linked to the character evolution of the trait) using multidimensional optimization
void TraitRELAXManager::optimizeAlternativeCharacterModelMultiDim()
{
  ParameterList parametersToEstimate = traitRELAXLikelihoodFunction_->getCharacterLikelihoodFunction()->getModelForSite(0,0)->getParameters();

  double tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance", traitRELAXParameters_->getParams(), .000001);
  unsigned int nbEvalMax = ApplicationTools::getParameter<unsigned int>("optimization.max_number_f_eval", traitRELAXParameters_->getParams(), 1000000);

  string mhPath = ApplicationTools::getAFilePath("optimization.character.message_handler", traitRELAXParameters_->getParams(), false, false, "", false, "none", 1);
  OutputStream* messageHandler =
  (mhPath == "none") ? 0 :
  (mhPath == "std") ? ApplicationTools::message.get() :
  new StlOutputStream(new ofstream(mhPath.c_str(), ios::out));

  string prPath = ApplicationTools::getAFilePath("optimization.character.profiler", traitRELAXParameters_->getParams(), false, false, "", false, "none", 1);
  OutputStream* profiler =
  (prPath == "none") ? 0 :
  (prPath == "std") ? ApplicationTools::message.get() :
  new StlOutputStream(new ofstream(prPath.c_str(), ios::out));
  if (profiler)
    profiler->setPrecision(20);

  unsigned int optVerbose = ApplicationTools::getParameter<unsigned int>("optimization.verbose", traitRELAXParameters_->getParams(), 2, "", false, 1);

  string order  = ApplicationTools::getStringParameter("optimization.method.derivatives", traitRELAXParameters_->getParams(), "Gradient", "", true, 1);
  string optMethod;
  if (order == "Gradient")
  {
    optMethod = OptimizationTools::OPTIMIZATION_GRADIENT;
  }
  else if (order == "Newton")
  {
    optMethod = OptimizationTools::OPTIMIZATION_NEWTON;
  }
  else
    throw Exception("Option '" + order + "' is not known for 'optimization.method.derivatives'.");

  AbstractNumericalDerivative* fun = 0;

  // Build optimizer:
  Optimizer* optimizer = 0;
  if (optMethod == OptimizationTools::OPTIMIZATION_GRADIENT)
  {
    fun = new TwoPointsNumericalDerivative(traitRELAXLikelihoodFunction_);
    fun->setInterval(0.0000001);
    optimizer = new ConjugateGradientMultiDimensions(fun);
  }
  else if (optMethod == OptimizationTools::OPTIMIZATION_NEWTON)
  {
    fun = new ThreePointsNumericalDerivative(traitRELAXLikelihoodFunction_);
    fun->setInterval(0.0001);
    optimizer = new PseudoNewtonOptimizer(fun);
  }
  else
    throw Exception("OptimizationTools::optimizeBranchLengthsParameters. Unknown optimization method: " + optMethod);

  // Numerical derivatives:
  fun->setParametersToDerivate(parametersToEstimate.getParameterNames());

  optimizer->setVerbose(optVerbose);
  optimizer->setProfiler(profiler);
  optimizer->setMessageHandler(messageHandler);
  optimizer->setMaximumNumberOfEvaluations(nbEvalMax);
  optimizer->getStopCondition()->setTolerance(tolerance);

  // Optimize TreeLikelihood function:
  optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  optimizer->init(parametersToEstimate);
  optimizer->optimize();

  // We're done.
  unsigned int n = optimizer->getNumberOfEvaluations();
  ApplicationTools::displayResult("# TL evaluations", TextTools::toString(n));
  if (profiler) delete profiler;
  if (fun) delete fun; 
  if (optimizer) delete optimizer;
}

/******************************************************************************/

// optimizes he character model with respect to the joint likelihood function while leaving the sequence parameters fixed
void TraitRELAXManager::optimizeAlternativeCharacterModel()
{
    uint verbose = static_cast<unsigned int>(ApplicationTools::getIntParameter("verbose", traitRELAXParameters_->getParams(), 1));
    traitRELAXLikelihoodFunction_->setOptimizationScope(JointLikelihoodFunction::OptimizationScope(0)); // no optimization of sequece parameters is required at these stage
    string characterOptimizationMethod = ApplicationTools::getStringParameter("optimization.character_method", traitRELAXParameters_->getParams(), "FullD(derivatives=Newton)");
    map<string, double> optimalCharacterParameters;
    optimalCharacterParameters.clear();
    if (characterOptimizationMethod.compare("grid") == 0)
    {
      uint gridSize = static_cast<unsigned int>(ApplicationTools::getIntParameter("optimization.grid_size", traitRELAXParameters_->getParams(), 10));
      optimizeAlternativeCharacterModelByGrid(optimalCharacterParameters, verbose, gridSize);
    }
    else if (characterOptimizationMethod.compare("oneDimBrent") == 0)
    {
      optimizeAlternativeCharacterModelByOneDimBrent(optimalCharacterParameters, 1);
    }
    else if (characterOptimizationMethod.compare("powell") == 0)
    {
      optimizeAlternativeCharacterModelByPowell(optimalCharacterParameters, verbose);
    }
    else
    {
      optimizeAlternativeCharacterModelMultiDim();
    } 
}

/******************************************************************************/

// writes inference result to output files
void TraitRELAXManager::writeInferenceToOutput()
{
  // Write parameters to file:
  string parametersFile = ApplicationTools::getAFilePath("output.estimates", traitRELAXParameters_->getParams(), false, false, "none", 1);
  bool withAlias = ApplicationTools::getBooleanParameter("output.estimates.alias", traitRELAXParameters_->getParams(), true, "", true, 0);
  const SiteContainer*seqData = traitRELAXLikelihoodFunction_->getSequenceLikelihoodFunction()->getData();
  SubstitutionModelSet* seqModel = traitRELAXLikelihoodFunction_->getSequenceLikelihoodFunction()->getSubstitutionModelSet();
  TransitionModel* charModel = traitRELAXLikelihoodFunction_->getCharacterLikelihoodFunction()->getModel();

  ApplicationTools::displayResult("Output estimates to file", parametersFile);
  if (parametersFile != "none")
  {
    StlOutputStream out(new ofstream(parametersFile.c_str(), ios::out));
    out << "# Log likelihood = ";
    out.setPrecision(20) << (-traitRELAXLikelihoodFunction_->getValue());
    out.endLine();
    out << "# Number of coding sequence sites = ";
    out.setPrecision(20) << seqData->getNumberOfSites();
    out.endLine();
    out.endLine();
    out << "# Substitution model parameters:";
    out.endLine();
    PhylogeneticsApplicationTools::printParameters(charModel, out, 1, withAlias);
    out.endLine();
    PhylogeneticsApplicationTools::printParameters(seqModel, out, 1, withAlias);
    out.endLine();
  }

  // Write infos to file:
  string infosFile = ApplicationTools::getAFilePath("output.infos", traitRELAXParameters_->getParams(), false, false);
  if (infosFile != "none")
  {
    ApplicationTools::displayResult("Alignment information logfile", infosFile);
    ofstream out(infosFile.c_str(), ios::out);

    // Get the rate class with maximum posterior probability:
    vector<size_t> classes = traitRELAXLikelihoodFunction_->getSequenceLikelihoodFunction()->getRateClassWithMaxPostProbOfEachSite();

    // Get the posterior rate, i.e. rate averaged over all posterior probabilities:
    Vdouble rates = traitRELAXLikelihoodFunction_->getSequenceLikelihoodFunction()->getPosteriorRateOfEachSite();

    vector<string> colNames;
    colNames.push_back("Sequence.Sites");
    colNames.push_back("is.complete");
    colNames.push_back("is.constant");
    colNames.push_back("lnL");
    colNames.push_back("rc");
    colNames.push_back("pr");
    vector<string> row(6);
    DataTable* infos = new DataTable(colNames);
  
  vector<double> LoglBySite = traitRELAXLikelihoodFunction_->getLikelihoodForEachSite();
  
    for (unsigned int i = 0; i < seqData->getNumberOfSites(); i++)
    {
      double lnL = LoglBySite[i];
      const Site* currentSite = &seqData->getSite(i);
      int currentSitePosition = currentSite->getPosition();
      string isCompl = "NA";
      string isConst = "NA";
      try { isCompl = (SiteTools::isComplete(*currentSite) ? "1" : "0"); }
      catch(EmptySiteException& ex) {}
      try { isConst = (SiteTools::isConstant(*currentSite) ? "1" : "0"); }
      catch(EmptySiteException& ex) {}
      row[0] = (string("[" + TextTools::toString(currentSitePosition) + "]"));
      row[1] = isCompl;
      row[2] = isConst;
      row[3] = TextTools::toString(lnL);
      row[4] = TextTools::toString(classes[i]);
      row[5] = TextTools::toString(rates[i]);
      infos->addRow(row);
    }

    DataTable::write(*infos, out, "\t");

    delete infos;
  }
}

// fit the alternative model: sequencial optimization of the character model and then sequence model, given an expected history based on the character model
map<string, double> TraitRELAXManager::optimizeAlternativeModel()
{
    cout << "\n******************************* Alternative model fitting *******************************" << endl;
    unsigned int verbose = static_cast<unsigned int>(ApplicationTools::getIntParameter("optimization.verbose", traitRELAXParameters_->getParams(), 0));
    bool report = false;
    if (verbose)
      report = true;
  
    traitRELAXParameters_->startTimer();

    // switch settings to alternative model (-> link between the character and sequence models)
    traitRELAXLikelihoodFunction_->setHypothesis(JointLikelihoodFunction::Hypothesis(1));
    traitRELAXLikelihoodFunction_->setOptimizationScope(JointLikelihoodFunction::OptimizationScope(2));
    
    // starting point optimization: naximum parsimony based partition
    if (verbose)
      cout << "\nStarting point optimization: maximum parsimony partition" << endl;
    traitRELAXParameters_->startTimer();
    ParameterList emptyParametersList;
	  traitRELAXLikelihoodFunction_->fireParameterChanged(emptyParametersList);
	
    // two-step optimization procedure
    if (verbose)
      cout << "* Starting iterative two-layer optimzation of the alternative model *" << endl;
    
    traitRELAXLikelihoodFunction_->getSequenceLikelihoodFunction()->computeTreeLikelihood(); // compute the initial log likelihood
    traitRELAXLikelihoodFunction_->getModelParameters(report); // report final values and complete step
    traitRELAXParameters_->done();

    
    double currentAlternativeOverallLogL = -traitRELAXLikelihoodFunction_->getValue();
    double previousAlternativeOverallLogL = currentAlternativeOverallLogL;
    map<string, double> bestModelParameters = traitRELAXLikelihoodFunction_->getModelParameters(false);
    map<string, double> currentmodelParameters;
    double bestLogl = currentAlternativeOverallLogL;
    double optimizationCyclesNum = 0;
    double tolerance = 0.1;
    do
    {
      traitRELAXParameters_->startTimer();
      /* optimize the character model while fixing the sequence model with respect to the joint model */
      cout << "\n** Step 1: fix sequence model parameters, optimize character model parameters **\n" << endl;
      traitRELAXLikelihoodFunction_->setOptimizationScope(JointLikelihoodFunction::OptimizationScope(0)); // disable sequence model parameters optimzation for the first step
      optimizeAlternativeCharacterModel();

      /* report optimal character parameters */
      cout << "\n** Optimal character parameters: **\n" << endl;
      ApplicationTools::displayResult("Mu", TextTools::toString(traitRELAXLikelihoodFunction_->getParameterValue("TwoParameterBinary.mu")));
      ApplicationTools::displayResult("Pi0", TextTools::toString(traitRELAXLikelihoodFunction_->getParameterValue("TwoParameterBinary.pi0")));
      traitRELAXParameters_->done();

      /* optimize the sequence model while fixing the character model with respect to the joint model */
      cout << "\n** Step 2: fix character model parameters, optimize sequence model parameters **\n" << endl;
      traitRELAXParameters_->startTimer();
      traitRELAXLikelihoodFunction_->setOptimizationScope(JointLikelihoodFunction::OptimizationScope(2)); // enable sequence model parameters optimzation for the first step
      traitRELAXLikelihoodFunction_->fireParameterChanged(emptyParametersList);
      currentmodelParameters = traitRELAXLikelihoodFunction_->getModelParameters(false);
      previousAlternativeOverallLogL = currentAlternativeOverallLogL;
      currentAlternativeOverallLogL = -traitRELAXLikelihoodFunction_->getValue();
      traitRELAXParameters_->done();
      if (currentAlternativeOverallLogL > bestLogl)
      {
        bestModelParameters = currentmodelParameters;
        bestLogl = -bestModelParameters["Overall Log likelihood"];
      }
      ApplicationTools::displayResult("Optimization cycle", TextTools::toString(optimizationCyclesNum + 1));
      ApplicationTools::displayResult("Difference between current and previous joint log likelihood", TextTools::toString(currentAlternativeOverallLogL - previousAlternativeOverallLogL));
      optimizationCyclesNum = optimizationCyclesNum + 1;
      cout << "\n\n************************************************************************************************\n\n" << endl;
    } while (currentAlternativeOverallLogL - previousAlternativeOverallLogL > tolerance);

  // report the optimal log likelihood and parameters of the alternative model
  for (map<string, double>::iterator it = bestModelParameters.begin(); it != bestModelParameters.end(); it++)
  {
    if ((it->first.find("Log likelihood") == std::string::npos))
    {
      ApplicationTools::displayResult(it->first, TextTools::toString(it->second, 15));
    }
  }
  reportParameters(bestModelParameters);
  traitRELAXParameters_->startTimer();
  ApplicationTools::displayResult("Number of optimization cycles", TextTools::toString(optimizationCyclesNum));
  alternativeLogl_ = -traitRELAXLikelihoodFunction_->getValue();
  writeInferenceToOutput();

  return bestModelParameters;
}

/******************************************************************************/
/***************************** Statistical testing ****************************/
/******************************************************************************/

// perfroms a traditional Llikelihood ratio test according to a chi2 null distribution of likelihood ratios with 1 degree of freedom
bool TraitRELAXManager::LRT()
{
    double LRTStatistic = 2*(alternativeLogl_ - nullLogl_);
    if (LRTStatistic > 3.84) // this fixed value corresponds to the threshold of a chi2 distribution with one degree of freedom
    {
      return true;
    }
    return false;
}

//simulates sequence data under the null hypothesis given the existing character data and the parameters obtained from the null optimization
//these data can then be analyzed using TraitRELAX to obtain a distribution of a null LL distribution
void TraitRELAXManager::ParametricBoostrap(map <string, double> simulationParameters, size_t numOfReplicates)
{
  // extract character parameters
  traitRELAXParameters_->getParam("character_model.mu") = TextTools::toString(simulationParameters["TwoParameterBinary.mu"]);
  traitRELAXParameters_->getParam("character_model.pi0") = TextTools::toString(simulationParameters["TwoParameterBinary.pi0"]);

  // extract sequence model parameters
  traitRELAXParameters_->getParam("sequence_model.kappa") = TextTools::toString(simulationParameters["RELAX.kappa_1"]);
  traitRELAXParameters_->getParam("sequence_model.omega0") = TextTools::toString(simulationParameters["RELAX.p_1"] * simulationParameters["RELAX.omega1_1"]);
  traitRELAXParameters_->getParam("sequence_model.p0") = TextTools::toString(simulationParameters["RELAX.theta1_1"]);
  traitRELAXParameters_->getParam("sequence_model.omega1") = TextTools::toString(simulationParameters["RELAX.omega1_1"]);
  traitRELAXParameters_->getParam("sequence_model.p1") = TextTools::toString((1-simulationParameters["RELAX.theta1_1"]) * simulationParameters["RELAX.theta2_1"]);
  traitRELAXParameters_->getParam("sequence_model.omega2") = TextTools::toString(simulationParameters["RELAX.omega2_1"]);
  traitRELAXParameters_->getParam("sequence_model.k") = TextTools::toString("RELAX.k_1"); // should be 1 under the null hypothesis
  traitRELAXParameters_->getParam("sequence_model.frequencies") = "F3X4";
  traitRELAXParameters_->getParam("sequence_model.frequencies.codon_pos1") = TextTools::toString(simulationParameters["RELAX.1_Full.theta_1"]) + ", " + TextTools::toString(simulationParameters["RELAX.1_Full.theta1_1"]) + ", " + TextTools::toString(simulationParameters["RELAX.1_Full.theta2_1"]);
  traitRELAXParameters_->getParam("sequence_model.frequencies.codon_pos2") = TextTools::toString(simulationParameters["RELAX.2_Full.theta_1"]) + ", " + TextTools::toString(simulationParameters["RELAX.2_Full.theta1_1"]) + ", " + TextTools::toString(simulationParameters["RELAX.2_Full.theta2_1"]);
  traitRELAXParameters_->getParam("sequence_model.frequencies.codon_pos3") = TextTools::toString(simulationParameters["RELAX.3_Full.theta_1"]) + ", " + TextTools::toString(simulationParameters["RELAX.3_Full.theta1_1"]) + ", " + TextTools::toString(simulationParameters["RELAX.3_Full.theta2_1"]);

  // add to traitRELAXParameters_ the required fields for simulation
  traitRELAXParameters_->getParam("sequence.num_of_sites") = TextTools::toString(traitRELAXLikelihoodFunction_->getSequenceLikelihoodFunction()->getNumberOfSites());
  traitRELAXParameters_->getParam("output.sequence.format") = "Fasta";
  traitRELAXParameters_->getParam("output.tree.format") = "Newick";

  // simulate trait evolution under the constraint of the states at the leaves
  // in other words, create a single stochastic mapping and use it as the trait history
  HomogeneousTreeLikelihood* characterLikelihoodFunction = traitRELAXLikelihoodFunction_->getCharacterLikelihoodFunction();
  characterLikelihoodFunction->setParameterValue("TwoParameterBinary.mu", simulationParameters["TwoParameterBinary.mu"]);
  characterLikelihoodFunction->setParameterValue("TwoParameterBinary.pi0", simulationParameters["TwoParameterBinary.pi0"]);
  StochasticMapping* stocMapping = new StochasticMapping(characterLikelihoodFunction, 1);
  vector<Tree*> mappings;
  stocMapping->generateStochasticMapping(mappings);
  Tree* traitHistory = mappings[0];

  // write the simulated history
  string outputDir = ApplicationTools::getStringParameter("output.simulations", traitRELAXParameters_->getParams(), GetCurrentWorkingDir(), "", true, 1);
  writeMapping(traitHistory, outputDir + "/unlabeled_trait_history.nwk", outputDir + "/labeled_trait_history.nwk");

  // simulate sequence evolution based on the partition derived from the trait evolution (which is now available in simulationParams)
  for (size_t rep=0; rep<numOfReplicates; ++rep)
  {
    traitRELAXParameters_->getParam("sequence.data_path") = outputDir + "/sequence_data_" + TextTools::toString(rep) + ".fas";
    simulateSequenceEvolution(traitHistory);
  }

  // free
  delete traitHistory;
  delete stocMapping;
}

/******************************************************************************/

// performs a statistical test according to the choice of the user inserted into the parameters file
void TraitRELAXManager::test(map<string, double> nullParameters, map<string, double> alternativeParameters)
{
  cout << "\n****************************** Performing statistical test ******************************" << endl;
 
  string testMethod = ApplicationTools::getStringParameter("statisical.test", traitRELAXParameters_->getParams(), "traditional LRT", "", true, 1);
  if (testMethod.compare("traditional LRT") == 0)
  {
    cout << "Tranditional Likelihood ratio test was selected..." << endl;
    bool isRejected = LRT();
    if (isRejected)
      cout << "Null hypothesis rejected with selection intensity k = " << alternativeParameters["RELAX.k_2"] << endl;
    else
      cout << "Null hypothesis not rejected " << endl;
  }
  else
  {
    cout << "Simulating data under the null hypothesis for parametric bootstrapping..." << endl;
    size_t numOfReplicates = ApplicationTools::getIntParameter("simulations.number_of_replicates", traitRELAXParameters_->getParams(), 100, "", true, 1);
    ParametricBoostrap(nullParameters, numOfReplicates);
  }
}

/* ################################################################## */
/* ####################### Simulation functions ##################### */
/* ################################################################## */

//removes substring from a string
void TraitRELAXManager::eraseSubStr(std::string & mainStr, const std::string & toErase)
{
    size_t pos = std::string::npos;
    // Search for the substring in string in a loop untill nothing is found
    while ((pos  = mainStr.find(toErase) )!= std::string::npos)
    {
        // If found then erase it from string
        mainStr.erase(pos, toErase.length());
    }
}

/******************************************************************************/

// get current working directory
string TraitRELAXManager::GetCurrentWorkingDir() 
{
  char buff[FILENAME_MAX];
  getcwd( buff, FILENAME_MAX );
  string current_working_dir(buff);
  return current_working_dir;
}

/* ################################################################## */
/* #################### Trait data simulator ######################## */
/* ################################################################## */

// since the function is fitted to true history, the mutation path may exceed the original branch length, in which case the dutration until the last event should be reduced to fit the original branch length
void TraitRELAXManager::updateBranchMapping(Node* son, const MutationPath& branchMapping, size_t initial_count)
{ 
    static size_t nodesCounter_ = initial_count;
    const vector<size_t> states = branchMapping.getStates();
    const VDouble times = branchMapping.getTimes();
    Node* curNode = son; // builds the branch history bottom to top (from last event to first event)
    double origBranchLength = son->getDistanceToFather();
    Node* nextNode;
    int eventsNum = static_cast<int>(branchMapping.getNumberOfEvents());
    if (eventsNum == 0)  // if there are no events -> return nothing
    {
        return;
    }
    else
    {
        // update the branch length of the node to be as per the time of the last event
        double duration = 0;
        for (int i=eventsNum-1; i>-1; --i) 
        { // add a new node to represent the transition
            nodesCounter_ += 1;
            const string name = "_mappingInternal" + TextTools::toString(nodesCounter_) + "_";
            nextNode = new Node(static_cast<int>(nodesCounter_), name); // nodes counter isn't incremented properly
            if (states[i] == 0)
            {
                StochasticMapping::setNodeState(nextNode,1);
            }
            else 
            {
                StochasticMapping::setNodeState(nextNode,0);
            }
            if (i == 0)
            {
                duration = times[i];    
            }
            else
            {
                duration = times[i]-times[i-1]; // the duration is the gap between the two adjacent times
            }
            nextNode->setDistanceToFather(duration); // according to the simulator, the documented time is the accumulated time until the transition and NOT the time since the last transition
            
            // set the father to no longer be the father of curNode
            if (curNode->hasFather()) // as long as we haven't reached the root -> upate the father of the current node
            {
                Node* originalFather = curNode->getFather();
                originalFather->removeSon(curNode); // also removes originalFather as the father of curNode
                // set nextNode to be the new father of curNode
                curNode->setFather(nextNode); // also adds curNode to the sons of nextNode
                // set curNode's original father ot be the father of nextNode
                nextNode->setFather(originalFather); // also adds nextNode to the sons of originalFather - or not? make sure this is the father at all times
                curNode = nextNode;
            }
        }
        // if the sum of distances doesn't add up to the length of the original branch, add an event or change the duration of the last event, depending on the last assigned state
        double lastEventDuration = origBranchLength - times[eventsNum-1]; 
        if (lastEventDuration > 0)
        {
            son->setDistanceToFather(lastEventDuration);
        }
    }
    return;
}

/*********************************************************************/

void TraitRELAXManager::addStatesToNodesNames(Tree* mapping)
{
    vector<Node*> nodes = (dynamic_cast<TreeTemplate<Node>*>(mapping))->getNodes(); 
    for (int i=0; i < static_cast<int>(nodes.size()); i++) 
    {
        string name = nodes[i]->getName();
        size_t state = StochasticMapping::getNodeState(nodes[i]);
        nodes[i]->setName(name + "{" + TextTools::toString(state) + "}");
    }
}

/*********************************************************************/

void TraitRELAXManager::removeStatesFromNodesNames(Tree* mapping)
{
    vector<Node*> nodes = (dynamic_cast<TreeTemplate<Node>*>(mapping))->getNodes(); 
    for (int i=0; i < static_cast<int>(nodes.size()); i++) 
    {
        string name = nodes[i]->getName();
        eraseSubStr(name, "{0}");
        eraseSubStr(name, "{1}");
        nodes[i]->setName(name);
    }
}

/*********************************************************************/

// translate the simulation data into a tree format, similar to the one given by the stochastic mappings
// in the results, there is a mutations paths field that can be returned for each node getMutationPath(int nodeId)
// then, it is enough to use the function I built in StochasticMapping
Tree* TraitRELAXManager::extractMapping(RASiteSimulationResult* simulationData)
{
    Tree* history = tree_->clone();
    giveNamesToInternalNodes(history);
    vector<Node*> nodes = (dynamic_cast<TreeTemplate<Node>*>(history))->getNodes();
    vector<Node*> debugNodes;
    size_t state;
    // set the ancestral states for all the nodes
    for (size_t j=0; j<nodes.size(); ++j)
    {
        // bypass bug in DetailedSiteSimulator.h line 127 - for root ancestral state, need to access ancestralStates[0] but can't because its index in indexes_ is 0 -> so line 127 end up accessing index 1
        // to bypass, call function in line 125 that accesses the argument index directly
        string nodeName = nodes[j]->getName();
        if (nodes[j]->getId() == history->getRootId())
        {
            state = simulationData->getRootAncestralState();
        }
        else
        {
            state = simulationData->getAncestralState(nodes[j]->getId());
        }
        StochasticMapping::setNodeState(nodes[j], state);
    }
    for (size_t i=0; i<nodes.size(); ++i)
    {
        string name = nodes[i]->getName();
        if (nodes[i]->hasFather())
        {
            // extract the path that ends in the visited node and update it on the tree
            const MutationPath branchHistory = simulationData->getMutationPath(nodes[i]->getId());
            updateBranchMapping(nodes[i], branchHistory, nodes.size()-1);
        } 
    }
    return history;
}

/* ################################################################## */

// the function verifies that the simulated trait evolution stricktly corresponds to the input phylgoeny strcuture
void TraitRELAXManager::checkMapping(Tree* mapping)
{
    // make sure both trees have the same size
    VDouble tree_BranchLengths = dynamic_cast<TreeTemplate<Node>*>(tree_)->getBranchLengths();
    double tree_Size = 0;
    for (size_t b=0; b<tree_BranchLengths.size(); b++)
    {
        tree_Size = tree_Size + tree_BranchLengths[b];
    }

    VDouble traitHistoryBranchLengths = dynamic_cast<TreeTemplate<Node>*>(mapping)->getBranchLengths();
    double traitHistorySize = 0;
    for (size_t b=0; b < traitHistoryBranchLengths.size(); b++)
    {
        traitHistorySize = traitHistorySize + traitHistoryBranchLengths[b];
    }

    if (abs(tree_Size - traitHistorySize) > 0.01)
    {
        throw Exception("true history has different tbl from base tree");
    }

    string nodeName;
    string mappingNodeName = "mapping";
    vector<Node*> traitHistoryNodes = (dynamic_cast<TreeTemplate<Node>*>(mapping))->getNodes();
    
    // make sure that each son of a base node has the same state as the base node
    Node* son;
    size_t nodeState, sonState;
    for (size_t n=0; n<traitHistoryNodes.size(); n++)
    {
        nodeName = traitHistoryNodes[n]->getName();
        // check if the node is a base node (either baseInternal of leaf)
        if (nodeName.find(mappingNodeName) == std::string::npos)
        {
            nodeState = StochasticMapping::getNodeState(traitHistoryNodes[n]);
            size_t numberOfSons = traitHistoryNodes[n]->getNumberOfSons();
            for (size_t s=0; s<numberOfSons; s++)
            {
                son = traitHistoryNodes[n]->getSon(s);
                sonState = StochasticMapping::getNodeState(son);
                if (sonState != nodeState)
                {
                    // if the son is a mapping node, its state should be the same as the state of its father base node
                    if (son->getName().find(mappingNodeName) != std::string::npos)
                    {
                        throw Exception("son state doesn't match parent state at node: " + son->getName());
                    }
                }
            }
        }  
    }

    // make sure that each parent of a base node has the same state as the base node
    size_t fatherState;
    Node* father;
    for(size_t i=0; i<traitHistoryNodes.size(); i++)
    {
        nodeName = traitHistoryNodes[i]->getName();
        // check if the node is a base node (either baseInternal of leaf)
        if (traitHistoryNodes[i]->getId() != mapping->getRootId())
        {
            if (!(nodeName.find(mappingNodeName) != std::string::npos))
            {
                nodeState = StochasticMapping::getNodeState(traitHistoryNodes[i]);
                father = traitHistoryNodes[i]-> getFather();
                fatherState = StochasticMapping::getNodeState(father);
                if (fatherState == nodeState)
                {
                    // if the parent is a mappimg node, its state must be differnt from the state of its base tree son
                    if (father->getName().find(mappingNodeName) != std::string::npos)
                    {
                        throw Exception("father state matches parent state at node: " + father->getName());
                    }
                    
                }
            }
        }  
    }
}

/* ################################################################## */

void TraitRELAXManager::derivePartition(Tree* mapping)
{
    vector<const Node*> nodes = (dynamic_cast<const TreeTemplate<Node>*>(mapping))->getNodes();
    // set assignment to branches
    map<size_t,string> stateToNodes;
    stateToNodes[0] = "";
    stateToNodes[1] = "";
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        if (!mapping->isRoot(nodes[i]->getId()))
            stateToNodes[StochasticMapping::getNodeState(nodes[i])] += TextTools::toString(nodes[i]->getId()) + ",";
    }
	// note: these labels will not be conserved upon reading of the trait evolutionary tree from scratch as node IDs will be provdied in post-order rather than added to a base tree like done during simulation
  traitRELAXParameters_->getParam("model1.nodes_id") = stateToNodes[0].substr(0, stateToNodes[0].length()-1);
  traitRELAXParameters_->getParam("model2.nodes_id") = stateToNodes[1].substr(0, stateToNodes[1].length()-1);
}

/* ################################################################## */

// write the mapping to two files (the first with the structure only, the seocnf with states labeling as well) and add the nodes mapping to models to the parameters instance accordingly
void TraitRELAXManager::writeMapping(Tree* mapping, string unlabeled_mapping_path, string labeled_mapping_path)
{
  std::map<std::string, std::string> writeParams;
  
  // write the true history to a file before deleting it
  writeParams["output.tree.file"]  = unlabeled_mapping_path;
  writeParams["output.tree.format"] = ApplicationTools::getStringParameter("output.tree.format", traitRELAXParameters_->getParams(), "Newick", "", true, 1);
  
  // write the tree without labels
  PhylogeneticsApplicationTools::writeTree(*mapping, writeParams);
  
  // now, add the label of each internal node in the updated history to the node's name
  addStatesToNodesNames(mapping);
  
  // write the tree with labels
  writeParams["output.tree.file"]  = labeled_mapping_path;
  PhylogeneticsApplicationTools::writeTree(*mapping, writeParams);
  
  // derive the partition of the trait history tree and write it to traitRELAXParameters_
  derivePartition(mapping);
  
  // now, remove the labels from the node names
  removeStatesFromNodesNames(mapping);
}


// simulate trait evolution along the provided phylogeny and write to output files: character data, simulated history in a newick format, and write to the BppApplication instance a map of nodes to model1 (BG) and model2 (FG)
Tree* TraitRELAXManager::simulateTraitEvolution()
{
    // extract from the parameters files the simulation parameters (default parameter values are set according to the manuscript)
    double mu = ApplicationTools::getDoubleParameter("character_model.mu", traitRELAXParameters_->getParams(), 8, "", true, 2);
    double pi0 = ApplicationTools::getDoubleParameter("character_model.pi0", traitRELAXParameters_->getParams(), 0.5, "", true, 2);

    // create a binary model
    if (!binaryAlphabet_)
		binaryAlphabet_ = new BinaryAlphabet();
    TwoParameterBinarySubstitutionModel* charModel = new TwoParameterBinarySubstitutionModel(binaryAlphabet_,mu, pi0); // second arguent stands for mu
    if (!rDist_)
		rDist_ = new ConstantRateDistribution();
    vector<string> seqNames = tree_->getLeavesNames();
    VectorSiteContainer charData(seqNames, binaryAlphabet_);
     
    // simulate character history using a simulator over a simple binary model
    NonHomogeneousSequenceSimulator* charSimulator = new NonHomogeneousSequenceSimulator(charModel, rDist_, tree_);

    RASiteSimulationResult* charResult = charSimulator->dSimulateSite();
    unique_ptr<Site> charSite(charResult->getSite(*charSimulator->getSubstitutionModelSet()->getModel(0)));
    charSite->setPosition(0);
    charData.addSite(*charSite, false);
    Tree* traitHistory = extractMapping(charResult);

    // check if the true history is legal
    checkMapping(traitHistory);

    // write the simulate character data
    map<string, string> writeParams;
    writeParams["output.sequence.file"] = ApplicationTools::getAFilePath("character.data_path", traitRELAXParameters_->getParams(), false, false, "", true, GetCurrentWorkingDir() + "/character_data.fas", 1);
    writeParams["output.sequence.format"] = ApplicationTools::getStringParameter("output.sequence.format", traitRELAXParameters_->getParams(), "Fasta", "", true, 1);
    SequenceApplicationTools::writeSequenceFile(*(dynamic_cast<SequenceContainer*>(&charData)), writeParams, "", false, 1);

    // write the true history to a file before deleting it
    string unlabled_history_path = ApplicationTools::getAFilePath("character.unlabeled_history_path", traitRELAXParameters_->getParams(), false, false, "", true, GetCurrentWorkingDir() + "/.unlabeled_trait_history.nwk", 1);
    string labeled_history_path  = ApplicationTools::getAFilePath("character.labeled_history_path", traitRELAXParameters_->getParams(), false, false, "", true, GetCurrentWorkingDir() + "/.labeled_trait_history.nwk", 1);
    writeMapping(traitHistory, unlabled_history_path, labeled_history_path);

    // free
    delete charModel;
    delete charSimulator;
    delete charResult;
    return traitHistory;
}

/* ################################################################## */
/* ################### Sequence data simulator ###################### */
/* ################################################################## */

//  sets the sequence model to simulate a partition of the sequence data with a specific omega
SubstitutionModelSet* TraitRELAXManager::setSequenceSubModel(const VectorSiteContainer* codon_data, double omega, double k)
{
  // process the other parameters to simulate with
  double kappa =  ApplicationTools::getDoubleParameter("sequence_model.kappa", traitRELAXParameters_->getParams(), 2, "", true, 0);
	string frequenciesModel = ApplicationTools::getStringParameter("sequence_model.frequencies", traitRELAXParameters_->getParams(), "F3X4", "", true, 0);
	string frequenciesInitalization = "";
	if (frequenciesModel.compare("F3X4") == 0)
	{
		vector<string> nucleotidesFrequenciesValues;
		vector<string> nucleotidesFrequenciesNames;
		string parName;
		nucleotidesFrequenciesNames.push_back("_Full.theta");
		nucleotidesFrequenciesNames.push_back("_Full.theta1");
		nucleotidesFrequenciesNames.push_back("_Full.theta2");
		for (size_t codonPos=1; codonPos<4; ++codonPos)
		{
			parName = "sequence_model.frequencies.codon_pos" + TextTools::toString(codonPos);

			nucleotidesFrequenciesValues = ApplicationTools::getVectorParameter<string>(parName, traitRELAXParameters_->getParams(), ',', "0.5,0.5,0.5", "", true, 0);
			for (size_t nucleotideIndex=0; nucleotideIndex<nucleotidesFrequenciesNames.size(); ++nucleotideIndex)
			{
				frequenciesInitalization += "," + TextTools::toString(codonPos) + nucleotidesFrequenciesNames[nucleotideIndex] + "=" + nucleotidesFrequenciesValues[nucleotideIndex];
			}
		}
	}
  traitRELAXParameters_->getParam("model1") = "YN98(kappa=" + TextTools::toString(kappa) + ",omega=" + TextTools::toString(omega) + ",frequencies=" + frequenciesModel + frequenciesInitalization + ")";
  // set the omega in model2 to be omega^k. The model is used for a specific partition the respective omegaWeight size and thus consists of a single omega parameter
  traitRELAXParameters_->getParam("model2") = "YN98(kappa=YN98.kappa_1,omega=" + TextTools::toString(pow(omega,k)) + ",frequencies=" + frequenciesModel + frequenciesInitalization + ")";
	traitRELAXParameters_->getParam("nonhomogeneous") = "general";
  traitRELAXParameters_->getParam("nonhomogeneous.number_of_models") = "2";
  traitRELAXParameters_->getParam("nonhomogeneous.stationarity") = "yes"; // constrain root frequencies to be the same as stationary (since RELAX is a time reversible model, this should not cause issues)

  // create the set of models
  if (!gCode_)
  {
    string codeDesc = ApplicationTools::getStringParameter("genetic_code", traitRELAXParameters_->getParams(), "Standard", "", true, true);
    gCode_ = SequenceApplicationTools::getGeneticCode(codonAlphabet_->getNucleicAlphabet(), codeDesc);
  }
  SubstitutionModelSet* modelSet = dynamic_cast<SubstitutionModelSet*>(PhylogeneticsApplicationTools::getSubstitutionModelSet(codonAlphabet_, gCode_, codon_data, traitRELAXParameters_->getParams()));
  return modelSet;
}

/* ################################################################## */

// adjusts the absoluate rates of the three independent models  to be as if they construct a single mixed model together
void TraitRELAXManager::homogenize(vector<SubstitutionModel*> models, vector<double>Probs)
{
    // step 0: compute synto and synfrom as - the minimal indices i,j of synonymous transition for which Q matrices of model 1 and model 2 have rate > 0 (that is, Q(i,j) > 0)
    // unforetunatly, the roginial data memebers are protected and thus must be comptated from scratch rather than extracted
    // emulated from YNGP_M2.cpp lines 110-121
    vector<int> supportedChars = models[0]->getAlphabetStates();
    size_t synfrom = 0;
	size_t synto = 0;
    for (synfrom = 1; synfrom < supportedChars.size(); ++synfrom)
    {
        for (synto = 0; synto < synfrom; ++synto)
        {
        if (gCode_->areSynonymous(supportedChars[synfrom], supportedChars[synto])
            && (models[0]->Qij(synfrom, synto) != 0)
            && (models[1]->Qij(synfrom, synto) != 0))
            break;
        }
        if (synto < synfrom)
        break;
    }
    if (synto == supportedChars.size())
        throw Exception("Impossible to find synonymous codons");

    // step 1: compute the vector of new rates as 1 / Q(synfrom, synto) for each model
    // emulated from YNGP_M2.cpp lines 136-142
    Vdouble Rates;
    Rates.push_back(1 / models[0]->Qij(synfrom, synto));
    Rates.push_back(1 / models[1]->Qij(synfrom, synto));
    Rates.push_back(1 / models[2]->Qij(synfrom, synto)); 

    // step 2: normalize the rates from step 1 as if all the models belong to the same mixture model: compute the weighted sum of rates
    // emulated from AbstractMixedSubstitutionModel.cpp lines 231-244
    double sum = 0;
    double sP = 0;
    for (unsigned int i = 0; i < Rates.size(); i++)
    {
        sum += Rates[i] * Probs[i];
        sP += Probs[i];
    }
    sum /= sP;

    // set the normalized rate of each model as its new rate
    for (unsigned int i=0; i<Rates.size(); i++)
    {
        Rates[i] *= 1 / sum;
        models[i]->setRate(Rates[i]);
    }
}

/* ################################################################## */

// simulate sequence evolution along the provided phylogeny, given the simulated trait evolution as partition, and write to output file the sequence data
void TraitRELAXManager::simulateSequenceEvolution(Tree* mapping)
{
    // extract simulation parameters
    size_t numOfSites = ApplicationTools::getIntParameter("sequence.num_of_sites", traitRELAXParameters_->getParams(), 300, "", true, 1);

    // sequence model parameters (the default parameter values are set according to the manuscript)
    vector<double>omegas;
    omegas.push_back(ApplicationTools::getDoubleParameter("sequence_model.omega0", traitRELAXParameters_->getParams(), 0.1, "", true, 1));
    omegas.push_back(ApplicationTools::getDoubleParameter("sequence_model.omega1", traitRELAXParameters_->getParams(), 0.8, "", true, 1));
    omegas.push_back(ApplicationTools::getDoubleParameter("sequence_model.omega2", traitRELAXParameters_->getParams(), 2, "", true, 1));
    vector<double> omegaWeights;
    omegaWeights.push_back(ApplicationTools::getDoubleParameter("p0", traitRELAXParameters_->getParams(), 0.5, "", true, 1));
    omegaWeights.push_back(ApplicationTools::getDoubleParameter("p1", traitRELAXParameters_->getParams(), 0.4, "", true, 1));
    omegaWeights.push_back(1-omegaWeights[0]-omegaWeights[1]);
    // process k values to simulate with
    double k_value = ApplicationTools::getDoubleParameter("k",  traitRELAXParameters_->getParams(), 1, "", true, 1);
    
    
    // create auxiliary variables for model generation
    TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(mapping);
    if (!rDist_)
		rDist_ = new ConstantRateDistribution();
    vector<string> seqNames = mapping->getLeavesNames();

    // generate an empty dataset and use it to create RELAX
    if (!codonAlphabet_)
		setCodonAlphabet();
    VectorSiteContainer seqDataContainter(seqNames, codonAlphabet_);

    // simulate a codon alignment according to RELAX model
    // RELAX is a mixture of 3 non homogenous models: one for each omega. Thus, to simulate under RELAX, we choose for each site the model to simulate it under, based on multimodel distrubution over the 3 models
    // create the 3 branch models: one per omega (site category)
    SubstitutionModelSet* seqModel_1 = setSequenceSubModel(&seqDataContainter, omegas[0], k_value);
    SubstitutionModelSet* seqModel_2 = setSequenceSubModel(&seqDataContainter, omegas[1], k_value);
    SubstitutionModelSet* seqModel_3 = setSequenceSubModel(&seqDataContainter, omegas[2], k_value);

    /* emulate homogenization */
    vector<SubstitutionModel*> bgModels;
    bgModels.push_back(dynamic_cast<SubstitutionModel*>(seqModel_1->getModel(0)));
    bgModels.push_back(dynamic_cast<SubstitutionModel*>(seqModel_2->getModel(0)));
    bgModels.push_back(dynamic_cast<SubstitutionModel*>(seqModel_3->getModel(0)));
    homogenize(bgModels, omegaWeights);
    vector<SubstitutionModel*> fgModels;
    fgModels.push_back(dynamic_cast<SubstitutionModel*>(seqModel_1->getModel(1)));
    fgModels.push_back(dynamic_cast<SubstitutionModel*>(seqModel_2->getModel(1)));
    fgModels.push_back(dynamic_cast<SubstitutionModel*>(seqModel_3->getModel(1)));
    homogenize(fgModels, omegaWeights);

    // simulate sites
    NonHomogeneousSequenceSimulator* seqSimulator_1 = new NonHomogeneousSequenceSimulator(seqModel_1, rDist_, ttree);
    NonHomogeneousSequenceSimulator* seqSimulator_2 = new NonHomogeneousSequenceSimulator(seqModel_2, rDist_, ttree);
    NonHomogeneousSequenceSimulator* seqSimulator_3 = new NonHomogeneousSequenceSimulator(seqModel_3, rDist_, ttree);
    VectorSiteContainer* seqData = new VectorSiteContainer(seqNames, codonAlphabet_);
    vector<size_t> omegaOrders;   // holds 1,2,3

    // for debugging purposes - track how many sites were simulated each omega category
    map <size_t, size_t> categroyToNumOfSites;
    for (size_t i=1; i<4; ++i)
    {
        categroyToNumOfSites[i] = 0;
    }

    for (size_t site=0; site<numOfSites; ++site)
    {

    RASiteSimulationResult* seqResult;
    Site* seqSite; 
    
    omegaOrders.clear();
    omegaOrders.push_back(1);
    omegaOrders.push_back(2);
    omegaOrders.push_back(3);
    
    size_t modelOfSite = RandomTools::pickOne(omegaOrders, omegaWeights, true);
    categroyToNumOfSites[modelOfSite] ++;  // for debugging purposes
    if (modelOfSite == 1)
    {
        seqResult = seqSimulator_1->dSimulateSite();
        seqSite = seqResult->getSite(*seqSimulator_1->getSubstitutionModelSet()->getModel(0));
    }
    else if (modelOfSite == 2)
    {
        seqResult = seqSimulator_2->dSimulateSite();
        seqSite = seqResult->getSite(*seqSimulator_2->getSubstitutionModelSet()->getModel(0));
    }
    else
    {
        seqResult = seqSimulator_3->dSimulateSite();
        seqSite = seqResult->getSite(*seqSimulator_3->getSubstitutionModelSet()->getModel(0));
    }
    seqSite->setPosition(static_cast<int>(site));
    seqData->addSite(*seqSite, false);
    delete seqSite;
    delete seqResult;
    }

    // for debugging purposes
    for (size_t i=1; i<4; ++i)
    {
        cout << "number of sites simulated under category " << i << " = " << categroyToNumOfSites[i] << endl;
    }
    
    // write the simulated sequence data
    map<string, string> writeParams;
    writeParams["output.sequence.file"] = ApplicationTools::getAFilePath("sequence.data_path", traitRELAXParameters_->getParams(), false, false, "", true, GetCurrentWorkingDir() + "/.sequence_data.fas.nwk", 1);
    writeParams["output.sequence.format"] = ApplicationTools::getStringParameter("output.sequence.format", traitRELAXParameters_->getParams(), "Fasta", "", true, 1);
    SequenceApplicationTools::writeSequenceFile(*(dynamic_cast<SequenceContainer*>(seqData)), writeParams, "", false, 1);

    // free
    delete seqData;
    delete seqModel_1;
    delete seqModel_2;
    delete seqModel_3;
    delete seqSimulator_1;
    delete seqSimulator_2;
    delete seqSimulator_3;
}

/* ################################################################## */

// simulates character and sequence data under the TraitELAX model
void TraitRELAXManager::simulate()
{
    // process the base tree
    setTree();

    // simulate trait evolution
    Tree* traitHistory = simulateTraitEvolution();

    // simulate sequence evolution based on the partition derived from the trait evolution (which is now available in simulationParams)
    simulateSequenceEvolution(traitHistory);

	// free
    delete traitHistory;
}

/* ################################################################## */