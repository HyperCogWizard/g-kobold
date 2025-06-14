#ifndef CELLULAR_VALIDATION_H
#define CELLULAR_VALIDATION_H

#include "../core/atomspace.h"
#include "../core/psystem.h"
#include <vector>
#include <map>
#include <string>

namespace cellular {

/**
 * Experimental data types for validation
 */
enum class ExperimentalDataType {
    CONCENTRATION_TIME_SERIES,    // Concentration vs time
    FLUX_MEASUREMENTS,           // Metabolic flux data
    ENZYME_KINETICS,             // Enzyme activity measurements
    MEMBRANE_POTENTIAL,          // Electrophysiology data
    PROTEIN_INTERACTIONS,        // Protein-protein interaction data
    GENE_EXPRESSION,             // mRNA/protein expression levels
    METABOLOMICS,                // Comprehensive metabolite profiling
    PROTEOMICS                   // Protein abundance data
};

/**
 * Cell biology experimental database
 */
class CellBiologyDatabase {
public:
    CellBiologyDatabase();
    
    // Standard cellular parameters from literature
    void loadStandardParameters();
    
    // Glucose metabolism (glycolysis, gluconeogenesis)
    std::map<std::string, double> getGlycosisKinetics() const;
    std::map<std::string, double> getTCAcycleKinetics() const;
    std::map<std::string, double> getElectronTransportKinetics() const;
    
    // Membrane transport
    std::map<std::string, double> getGlucoseTransporterKinetics() const;
    std::map<std::string, double> getSodiumPotassiumPumpKinetics() const;
    std::map<std::string, double> getCalciumChannelKinetics() const;
    
    // Signaling pathways
    std::map<std::string, double> getInsulinSignalingKinetics() const;
    std::map<std::string, double> getGPCRSignalingKinetics() const;
    std::map<std::string, double> getMapKinaseKinetics() const;
    
    // Cellular concentrations (typical values)
    std::map<std::string, double> getTypicalCytoplasmicConcentrations() const;
    std::map<std::string, double> getTypicalMitochondrialConcentrations() const;
    std::map<std::string, double> getTypicalNuclearConcentrations() const;
    
    // Organelle properties
    std::map<std::string, double> getMitochondrialProperties() const;
    std::map<std::string, double> getERProperties() const;
    std::map<std::string, double> getGolgiProperties() const;
    
    // Cell type specific parameters
    std::map<std::string, double> getHepatocyteParameters() const;
    std::map<std::string, double> getNeuronParameters() const;
    std::map<std::string, double> getCardiomyocyteParameters() const;
    
private:
    std::map<std::string, std::map<std::string, double>> experimental_parameters_;
    
    void loadGlycolysisData();
    void loadRespirationData();
    void loadSignalingData();
    void loadMembraneData();
    void loadOrganelleData();
};

/**
 * Model validation against experimental cell biology
 */
class ExperimentalValidator {
public:
    ExperimentalValidator();
    
    // Add experimental datasets
    void addExperimentalDataset(const std::string& name, 
                               ExperimentalDataType type,
                               const std::vector<double>& data,
                               const std::map<std::string, double>& conditions);
    
    // Validation tests
    double validateGlucoseMetabolism(const AtomSpace& atomspace, const PSystem& psystem) const;
    double validateATPProduction(const AtomSpace& atomspace, const PSystem& psystem) const;
    double validateMembraneTransport(const AtomSpace& atomspace, const PSystem& psystem) const;
    double validateSignalingResponse(const AtomSpace& atomspace, const PSystem& psystem) const;
    double validateMitochondrialFunction(const AtomSpace& atomspace, const PSystem& psystem) const;
    
    // Comprehensive validation
    std::map<std::string, double> performComprehensiveValidation(
        const AtomSpace& atomspace, const PSystem& psystem) const;
    
    // Statistical analysis
    double calculateRSquared(const std::vector<double>& experimental,
                           const std::vector<double>& simulated) const;
    double calculateRMSE(const std::vector<double>& experimental,
                        const std::vector<double>& simulated) const;
    double calculateMAPE(const std::vector<double>& experimental,
                        const std::vector<double>& simulated) const;
    
    // Generate validation report
    void generateValidationReport(const std::string& filename,
                                 const std::map<std::string, double>& validation_results) const;
    
private:
    struct ExperimentalDataset {
        ExperimentalDataType type;
        std::vector<double> data;
        std::map<std::string, double> conditions;
        std::string description;
    };
    
    std::map<std::string, ExperimentalDataset> experimental_datasets_;
    CellBiologyDatabase database_;
    
    // Validation helper functions
    std::vector<double> extractSimulationData(const AtomSpace& atomspace,
                                             const PSystem& psystem,
                                             const std::string& variable) const;
    
    double compareWithLiterature(const std::string& parameter,
                                double simulated_value) const;
};

/**
 * Specific experimental validation tests
 */
class ValidationTests {
public:
    // Classic cell biology experiments
    static double validateWarburg_effect(const AtomSpace& atomspace, const PSystem& psystem);
    static double validatePasteur_effect(const AtomSpace& atomspace, const PSystem& psystem);
    static double validateCrabtree_effect(const AtomSpace& atomspace, const PSystem& psystem);
    
    // Metabolic validation
    static double validateGlucose_dose_response(const AtomSpace& atomspace, const PSystem& psystem,
                                               const std::vector<double>& glucose_concentrations);
    static double validateInsulin_signaling(const AtomSpace& atomspace, const PSystem& psystem,
                                           double insulin_concentration);
    static double validateOxygen_consumption(const AtomSpace& atomspace, const PSystem& psystem);
    
    // Membrane validation
    static double validateGoldman_Hodgkin_Katz(const AtomSpace& atomspace, const PSystem& psystem);
    static double validateNernst_potential(const AtomSpace& atomspace, const PSystem& psystem,
                                          const std::string& ion);
    
    // Kinetic validation
    static double validateMichaelis_Menten_kinetics(const AtomSpace& atomspace, const PSystem& psystem,
                                                   const std::string& enzyme);
    static double validateHill_cooperativity(const AtomSpace& atomspace, const PSystem& psystem,
                                            const std::string& protein);
    
    // Organelle validation
    static double validateMitochondrial_respiration(const AtomSpace& atomspace, const PSystem& psystem);
    static double validateER_calcium_dynamics(const AtomSpace& atomspace, const PSystem& psystem);
    static double validateGolgi_processing(const AtomSpace& atomspace, const PSystem& psystem);
};

/**
 * Benchmark datasets from cell biology literature
 */
class CellBiologyBenchmarks {
public:
    // Glucose metabolism benchmarks
    struct GlucoseMetabolismBenchmark {
        std::vector<double> glucose_concentrations = {0, 1, 2.5, 5, 10, 20}; // mM
        std::vector<double> glucose_uptake_rates = {0, 0.5, 1.0, 1.5, 1.8, 1.9}; // mM/min
        std::vector<double> lactate_production = {0, 0.1, 0.3, 0.5, 0.7, 0.8}; // mM/min
        double km_glucose = 5.0; // mM
        double vmax_glucose = 2.0; // mM/min
    };
    
    // Mitochondrial respiration benchmark
    struct RespirationBenchmark {
        double state_3_respiration = 100.0; // nmol O2/min/mg protein
        double state_4_respiration = 15.0;  // nmol O2/min/mg protein
        double respiratory_control_ratio = 6.7; // State 3 / State 4
        double p_o_ratio = 2.5; // ADP/O ratio
        std::vector<double> substrate_concentrations = {0, 1, 2, 5, 10}; // mM
        std::vector<double> respiration_rates = {15, 35, 55, 85, 100}; // nmol O2/min/mg
    };
    
    // Calcium signaling benchmark
    struct CalciumBenchmark {
        double resting_calcium = 0.1e-6; // M (100 nM)
        double peak_calcium = 1.0e-6;    // M (1 μM)
        double half_time_rise = 0.5;     // seconds
        double half_time_decay = 5.0;    // seconds
        std::vector<double> stimulus_strengths = {0, 1, 10, 100, 1000}; // nM agonist
        std::vector<double> calcium_responses = {0.1, 0.2, 0.5, 0.8, 1.0}; // μM peak
    };
    
    static GlucoseMetabolismBenchmark getGlucoseMetabolismBenchmark();
    static RespirationBenchmark getRespirationBenchmark();
    static CalciumBenchmark getCalciumBenchmark();
    
    // Generate benchmark test suite
    static std::map<std::string, std::vector<double>> generateBenchmarkSuite();
    
    // Validate simulation against all benchmarks
    static std::map<std::string, double> validateAgainstAllBenchmarks(
        const AtomSpace& atomspace, const PSystem& psystem);
};

} // namespace cellular

#endif // CELLULAR_VALIDATION_H