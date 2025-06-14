#ifndef CELLULAR_SIMULATION_SYSTEM_H
#define CELLULAR_SIMULATION_SYSTEM_H

#include "core/atomspace.h"
#include "core/psystem.h"
#include "membrane/membrane_dynamics.h"
#include "signaling/signaling_pathways.h"
#include "organelles/organelles.h"
#include "organelles/mitochondria.h"
#include "kinetics/kinetics_engine.h"
#include <memory>
#include <string>
#include <vector>
#include <map>

namespace cellular {

/**
 * Experimental validation data structure
 */
struct ExperimentalData {
    std::string experiment_name;
    std::string cell_type;
    std::vector<std::string> measured_variables;
    std::vector<std::vector<double>> time_series_data;
    std::vector<double> time_points;
    std::map<std::string, double> experimental_conditions;
    double confidence_interval = 0.95;
};

/**
 * Simulation parameters
 */
struct SimulationParameters {
    double total_simulation_time = 100.0;  // seconds
    double time_step = 0.01;                // seconds
    double output_interval = 1.0;          // seconds
    bool enable_stochastic = false;
    bool enable_fluctuations = false;
    bool enable_multiscale = false;
    std::string integration_method = "runge_kutta_4";
    double tolerance = 1e-6;
};

/**
 * Cell type specifications
 */
enum class CellType {
    GENERIC_HUMAN,
    HEPATOCYTE,
    NEURON,
    CARDIOMYOCYTE,
    FIBROBLAST,
    EPITHELIAL,
    IMMUNE_CELL,
    STEM_CELL
};

/**
 * Comprehensive human cell simulation system
 */
class CellularSimulationSystem {
public:
    CellularSimulationSystem(CellType cell_type = CellType::GENERIC_HUMAN);
    ~CellularSimulationSystem() = default;
    
    // System initialization
    void initializeCell(CellType cell_type);
    void setSimulationParameters(const SimulationParameters& params);
    
    // Component access
    AtomSpace& getAtomSpace() { return atomspace_; }
    PSystem& getPSystem() { return psystem_; }
    MembraneDynamics& getMembraneDynamics() { return membrane_dynamics_; }
    IntracellularSignaling& getSignaling() { return signaling_; }
    OrganelleSystem& getOrganelleSystem() { return organelle_system_; }
    KineticsEngine& getKineticsEngine() { return kinetics_engine_; }
    
    // Cell-specific setup
    void setupGenericHumanCell();
    void setupHepatocyte();
    void setupNeuron();
    void setupCardiomyocyte();
    
    // Environmental conditions
    void setExtracellularConditions(const std::map<std::string, double>& conditions);
    void setTemperature(double temperature) { temperature_ = temperature; }
    void setpH(double ph) { ph_ = ph; }
    void setOxygenLevel(double oxygen) { oxygen_level_ = oxygen; }
    
    // Stimulus and perturbations
    void addLigandStimulus(const std::string& ligand, double concentration, double start_time, double duration = -1);
    void addDrugTreatment(const std::string& drug, double concentration, double start_time, double duration = -1);
    void addGeneticPerturbation(const std::string& gene, double expression_change);
    void addMetabolicPerturbation(const std::string& metabolite, double concentration_change);
    
    // Simulation control
    void runSimulation();
    void runTimeCourse(double total_time, double output_interval = 1.0);
    void runSteadyStateAnalysis();
    void runParameterScan(const std::string& parameter, double min_val, double max_val, int steps);
    
    // Analysis and output
    std::map<std::string, std::vector<double>> getTimeCourseData() const;
    std::map<std::string, double> getCurrentState() const;
    std::map<std::string, double> getMetabolicFluxes() const;
    std::map<std::string, double> getSignalingActivities() const;
    
    // Validation
    void addExperimentalData(const ExperimentalData& data);
    double validateAgainstExperiment(const std::string& experiment_name) const;
    std::map<std::string, double> calculateModelFitMetrics() const;
    
    // Outputs and visualization
    void exportResults(const std::string& filename, const std::string& format = "csv") const;
    void generateReport(const std::string& filename) const;
    
    // Mathematical framework export
    void exportMathematicalFramework(const std::string& filename) const;
    void exportPSystemRules(const std::string& filename) const;
    void exportAtomSpaceStructure(const std::string& filename) const;
    
    // Model introspection
    std::vector<std::string> listAllMolecules() const;
    std::vector<std::string> listAllReactions() const;
    std::vector<std::string> listAllPathways() const;
    std::map<std::string, double> getCompartmentVolumes() const;
    
    // Debug and diagnostics
    void enableLogging(bool enable) { logging_enabled_ = enable; }
    void setLogLevel(const std::string& level) { log_level_ = level; }
    std::string getDiagnosticInfo() const;
    
private:
    // Core components
    AtomSpace atomspace_;
    PSystem psystem_;
    MembraneDynamics membrane_dynamics_;
    IntracellularSignaling signaling_;
    OrganelleSystem organelle_system_;
    KineticsEngine kinetics_engine_;
    
    // Cell properties
    CellType cell_type_;
    double temperature_ = 310.15;  // 37°C in Kelvin
    double ph_ = 7.4;
    double oxygen_level_ = 0.21;   // 21% oxygen
    
    // Simulation state
    SimulationParameters sim_params_;
    double current_time_ = 0.0;
    std::map<std::string, std::vector<double>> time_course_data_;
    std::vector<double> time_points_;
    
    // Environmental and stimulus
    std::map<std::string, double> extracellular_conditions_;
    
    struct StimulusEvent {
        std::string type;
        std::string target;
        double concentration;
        double start_time;
        double duration;
        bool active;
    };
    std::vector<StimulusEvent> stimulus_events_;
    
    // Experimental validation
    std::vector<ExperimentalData> experimental_data_;
    
    // Logging and diagnostics
    bool logging_enabled_ = false;
    std::string log_level_ = "INFO";
    
    // Initialization helpers
    void initializeBasicMetabolism();
    void initializeSignalingPathways();
    void initializeTransportSystems();
    void initializeOrganelles();
    
    // Cell-type specific initialization
    void addHepatocyteSpecificFeatures();
    void addNeuronSpecificFeatures();
    void addCardiomyocyteSpecificFeatures();
    
    // Simulation helpers
    void updateStimuli(double current_time);
    void recordTimePoint();
    void logMessage(const std::string& level, const std::string& message) const;
    
    // Analysis helpers
    double calculateCorrelation(const std::vector<double>& experimental, 
                               const std::vector<double>& simulated) const;
    double calculateRMSE(const std::vector<double>& experimental, 
                        const std::vector<double>& simulated) const;
};

/**
 * Factory for creating preconfigured cellular simulation systems
 */
class CellularSimulationFactory {
public:
    static std::unique_ptr<CellularSimulationSystem> createHumanHepatocyte();
    static std::unique_ptr<CellularSimulationSystem> createHumanNeuron();
    static std::unique_ptr<CellularSimulationSystem> createHumanCardiomyocyte();
    static std::unique_ptr<CellularSimulationSystem> createGenericHumanCell();
    
    // Disease models
    static std::unique_ptr<CellularSimulationSystem> createDiabeticHepatocyte();
    static std::unique_ptr<CellularSimulationSystem> createAlzheimersNeuron();
    static std::unique_ptr<CellularSimulationSystem> createIschemiicCardiomyocyte();
    
    // Cancer models
    static std::unique_ptr<CellularSimulationSystem> createCancerCell(const std::string& cancer_type);
    
private:
    static void configureCellTypeSpecifics(CellularSimulationSystem& system, CellType type);
    static void addDiseaseModifications(CellularSimulationSystem& system, const std::string& disease);
};

/**
 * Batch simulation and analysis system
 */
class BatchSimulationSystem {
public:
    BatchSimulationSystem();
    
    // Batch operations
    void addSimulation(std::unique_ptr<CellularSimulationSystem> system, const std::string& name);
    void runAllSimulations();
    void runParameterSweep(const std::string& parameter, 
                          const std::vector<double>& parameter_values);
    
    // Comparative analysis
    std::map<std::string, std::map<std::string, double>> compareSimulations() const;
    void generateComparativeReport(const std::string& filename) const;
    
    // Statistical analysis
    void performSensitivityAnalysis(const std::vector<std::string>& parameters);
    void performUncertaintyQuantification(int monte_carlo_samples = 1000);
    
private:
    std::map<std::string, std::unique_ptr<CellularSimulationSystem>> simulations_;
    std::map<std::string, std::map<std::string, std::vector<double>>> batch_results_;
};

} // namespace cellular

#endif // CELLULAR_SIMULATION_SYSTEM_H