// From bpp-core:
#include <Bpp/Version.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Numeric/Function/BrentOneDimension.h>
#include <Bpp/Numeric/Function/PowellMultiDimensions.h>
#include <Bpp/Numeric/ParameterList.h>
#include <Bpp/Numeric/Function/BfgsMultiDimensions.h>
#include <Bpp/Numeric/Function/ReparametrizationFunctionWrapper.h>
#include <Bpp/Numeric/Function/ThreePointsNumericalDerivative.h>
#include <Bpp/Numeric/Function/ConjugateGradientMultiDimensions.h>
#include <Bpp/Numeric/Function/TwoPointsNumericalDerivative.h>
#include <Bpp/Numeric/Function/DownhillSimplexMethod.h>
#include <Bpp/Numeric/Random/RandomTools.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Io/BppOSequenceReaderFormat.h>
#include <Bpp/Seq/Io/BppOAlignmentReaderFormat.h>
#include <Bpp/Seq/Io/BppOSequenceWriterFormat.h>
#include <Bpp/Seq/Io/BppOAlignmentWriterFormat.h>
#include <Bpp/Seq/SequenceTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Node.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RNonHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RASTools.h>
#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Model/MixedTransitionModel.h>
#include <Bpp/Phyl/Model/TwoParameterBinarySubstitutionModel.h>
#include <Bpp/Phyl/Model/Protein/CoalaCore.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/FrequencySet/MvaFrequencySet.h>
#include <Bpp/Phyl/Model/FrequencySet/FrequencySet.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Io/BppOFrequencySetFormat.h>
#include <Bpp/Phyl/Mapping/StochasticMapping.h>
#include <Bpp/Phyl/Likelihood/JointLikelihoodFunction.h>
#include <Bpp/Phyl/Likelihood/PseudoNewtonOptimizer.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Simulation/DetailedSiteSimulator.h>
#include <Bpp/Phyl/Simulation/NonHomogeneousSequenceSimulator.h>
#include <Bpp/Phyl/Simulation/SequenceSimulationTools.h>

// From the STL:
#include <iostream>
#include <iomanip>
#include <limits>
#include <map>
#include <vector>
#include <stdio.h> 
#include <unistd.h>
 
using namespace bpp;
using namespace std;

typedef vector<vector<double>> VVDouble;
typedef vector<double> VDouble;
typedef unsigned int uint;

class TraitRELAXManager
{
    protected:
        BppApplication* traitRELAXParameters_;
        JointLikelihoodFunction*  traitRELAXLikelihoodFunction_;
        double nullLogl_;
        double alternativeLogl_;
        GeneticCode* gCode_;          // this is saved as a data member of that its pointer wouldn't be lost upon closure of init(), a stage in which it cannot be deleted yet
        //DiscreteDistribution* rDist_; // this is saved as a data member of that its pointer wouldn't be lost upon closure of init(), a stage in which it cannot be deleted yet

    public:

        // constructor
        TraitRELAXManager(BppApplication* parameters);
        
        // destructor
        ~TraitRELAXManager();

        // copy constructor
        TraitRELAXManager(const TraitRELAXManager& trm) : 
        traitRELAXParameters_(0), traitRELAXLikelihoodFunction_(), nullLogl_(0), alternativeLogl_(0), gCode_() //, rDist_()
        {   
            traitRELAXLikelihoodFunction_ = trm.traitRELAXLikelihoodFunction_; 
            traitRELAXParameters_ = trm.traitRELAXParameters_; 
            gCode_ = trm.gCode_;
            //rDist_ = trm.rDist_;
        }

        // assignment operator
        TraitRELAXManager& operator=(const TraitRELAXManager& trm)
        {
            traitRELAXParameters_ = trm.traitRELAXParameters_;
            traitRELAXLikelihoodFunction_ = trm.traitRELAXLikelihoodFunction_;
            gCode_ = trm.gCode_;
            //rDist_ = trm.rDist_;
            return *this;
        }

        // initialize the execution of traitrelax
        void init();

        // optimizes the joint model under the null hypothesis using independent optimizations of the characxter and sequence models
        // returns mapping of the model parameter names to their optimal value found during the search
        map<string, double> optimizeNullModel();

        // optimizes the joint model under the alternative hypothesis using a two step-optimzation procedure
        // first step is optimizing the character model parameters while leaving the sequence model parameters fixed
        // the second step is optimizing the sequence model parameters given a fixed set of character model parameters and accordingly a fixed trait history
        // returns mapping of the model parameter names to their optimal value found during the search
        map<string, double> optimizeAlternativeModel();

        // performs a statistical test according to the choice of the user inserted into the parameters file
        void test(map<string, double> nullParameters, map<string, double> alternativeParameters);

        // simulates character and sequence data under the TraitELAX model
        void simulate();    

    protected:

        // creates a codon alphabet
        const CodonAlphabet* getCodonAlphabet();

        // processes input character trait data
        VectorSiteContainer* processCharacterData(const BinaryAlphabet* alphabet);

        // processes input sequence data
        VectorSiteContainer* processSequenceData(const CodonAlphabet* codonAlphabet);

        // gives names to internal nodes for the trait histories internal nodes identification
        void giveNamesToInternalNodes(Tree* tree);

        // processes input tree
        Tree* processTree();

        // creates the character trait model
        TransitionModel* setCharacterModel(VectorSiteContainer* charData, const BinaryAlphabet* alphabet, Tree* tree, DRTreeParsimonyScore* mpData);

        // sets initial partition of the branches as input of the sequence model based on a maximum parsimony solution
        void setMpPartition(DRTreeParsimonyScore* mpData, const VectorSiteContainer* characterData, TransitionModel* characterModel, Tree* tree);

        // creates the sequence model
        MixedSubstitutionModelSet* setSequenceModel(const VectorSiteContainer* codon_data, const CodonAlphabet* codonAlphabet, DRTreeParsimonyScore* mpData, const VectorSiteContainer* characterData, TransitionModel* characterModel, Tree* tree);

        // reports parameter values and likelihood following optimization
        void reportParameters(map <string, double> parameterValues);
        
        // initalizes grid for character model optimization under the alternative hypothesis
        void initCharacterGrid(VVDouble &grid, double pi0Lb, double pi0Ub, double muLb, double muUb, uint gridSize);

        // obtains neighbors of a grid visited point
        VVDouble getNeighbors(VVDouble &grid, double currentMu, double currentPi0);

        // optimizes the character model under the alternative hypothesis (in which the sequence evolution is linked to the character evolution of the trait) using one dimensional grid optimization
        void optimizeAlternativeCharacterModelByGrid(map<string, double> &optimalValues, TransitionModel* characterModel, uint verbose = 1, uint gridSize = 10, bool firstCycle = true);

        // optimizes the character model under the alternative hypothesis (in which the sequence evolution is linked to the character evolution of the trait) using one dimensional brent optimization
        void optimizeAlternativeCharacterModelByOneDimBrent(map<string, double> &optimalValues, TransitionModel* characterModel, uint verbose = 1);

        // optimizes the character model under the alternative hypothesis (in which the sequence evolution is linked to the character evolution of the trait) using one dimensional powell optimization
        void optimizeAlternativeCharacterModelByPowell(map<string, double> &optimalValues, TransitionModel* characterModel, uint verbose = 1);

        // optimizes the character model under the alternative hypothesis (in which the sequence evolution is linked to the character evolution of the trait) using multidimensional optimization
        void optimizeAlternativeCharacterModelMultiDim();

        // optimizes the character model under the alternative hypothesis according to the selected method
        void optimizeAlternativeCharacterModel();

        // writes inference result to output files
        void writeInferenceToOutput();

        // perfroms a traditional Llikelihood ratio test according to a chi2 null distribution of likelihood ratios with 1 degree of freedom
        // returns true with null is rejected, else false
        bool LRT();

        //simulates sequence data under the null hypothesis given the existing character data and the parameters obtained from the null optimization
        //these data can then be analyzed using TraitRELAX to obtain a distribution of a null LL distribution
        void ParametricBoostrap(map<string, double> simulationParameters, size_t numderOfReplicates = 100);

        //removes substring from a string
        void eraseSubStr(std::string & mainStr, const std::string & toErase);

        // get curreent working directory
        string GetCurrentWorkingDir();

        // updates the mapping in a branch of a true history in the tree structure based on the given mutation process instance
        void updateBranchMapping(Node* son, const MutationPath& branchMapping, size_t initial_count);

        // adds to the nodes names their assigned states (upon writing the simulated history to a file). This will enable both visualization of the simulated histories and providing them as input partitionsd to hyphy
        void addStatesToNodesNames(Tree* mapping);

        // removes the states from the nodes names upon writing of the simulated history structure to a file. This will enable using the simulated hisotries as input to BIo++.
        void removeStatesFromNodesNames(Tree* mapping);

        // reconstructs the simulated history in a newick structure that can be given as input to the downstream sequence data simulator
        Tree* extractMapping(RASiteSimulationResult* simulationData, Tree* baseTree);

        // makes sure that the strcture of the simulated histories complies withh the input tree
        void checkMapping(Tree* mapping, Tree* baseTree);

        // extracts the partition induced by the simulated history
        void derivePartition(Tree* mapping);

        // write the mapping to two files (the first with the structure only, the seocnf with states labeling as well) and add the nodes mapping to models to the parameters instance accordingly
        void writeMapping(Tree* mapping, string unlabeled_mapping_path, string labeled_mapping_path);

        // simulates trait evolution along the provided phylogeny and write to output files: character data, simulated history in a newick format, and write to the BppApplication instance a map of nodes to model1 (BG) and model2 (FG)
        Tree* simulateTraitEvolution(Tree* tree);

        //  sets the sequence model to simulate a partition of the sequence data with a specific omega
        SubstitutionModelSet* setSequenceSubModel(const VectorSiteContainer* codon_data, const CodonAlphabet* codonAlphabet, double omega, double k);

        // adjusts the absoluate rates of the three independent models  to be as if they construct a single mixed model together
        void homogenize(vector<SubstitutionModel*> models, vector<double>Probs, const CodonAlphabet* codonAlphabet);

        // simulates sequence evolution along the provided phylogeny, given the simulated trait evolution as partition, and write to output file the sequence data
        void simulateSequenceEvolution(Tree* tree);
};
