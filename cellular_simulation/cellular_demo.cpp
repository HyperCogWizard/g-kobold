#include "cellular_simulation_system.h"
#include <iostream>
#include <iomanip>

using namespace cellular;

/**
 * Demonstration of the comprehensive cellular simulation system
 */
void demonstrate_cellular_simulation() {
    std::cout << "=== Cellular Simulation System Demonstration ===" << std::endl;
    std::cout << "Creating a generic human cell simulation..." << std::endl;
    
    // Create cellular simulation system
    auto cell_system = CellularSimulationFactory::createGenericHumanCell();
    
    // Set simulation parameters
    SimulationParameters params;
    params.total_simulation_time = 60.0;  // 60 seconds
    params.time_step = 0.1;               // 0.1 second steps
    params.output_interval = 1.0;         // Output every second
    params.enable_stochastic = false;
    params.integration_method = "runge_kutta_4";
    
    cell_system->setSimulationParameters(params);
    
    // Set physiological conditions
    std::map<std::string, double> extracellular_conditions = {
        {"glucose", 5.0},      // 5 mM glucose
        {"oxygen", 0.21},      // 21% oxygen
        {"Na+", 145.0},        // 145 mM sodium
        {"K+", 5.0},           // 5 mM potassium
        {"Ca2+", 2.0},         // 2 mM calcium
        {"Cl-", 110.0}         // 110 mM chloride
    };
    
    cell_system->setExtracellularConditions(extracellular_conditions);
    cell_system->setTemperature(310.15);  // 37°C
    cell_system->setpH(7.4);
    
    std::cout << "Initial cellular state established." << std::endl;
    std::cout << "Cell type: Generic Human Cell" << std::endl;
    std::cout << "Temperature: 37°C, pH: 7.4" << std::endl;
    std::cout << "Extracellular glucose: 5.0 mM" << std::endl;
    
    // Add a glucose stimulus
    std::cout << "\nAdding glucose stimulus at t=10s..." << std::endl;
    cell_system->addLigandStimulus("insulin", 1e-9, 10.0, 30.0);  // 1 nM insulin for 30s
    
    // Run simulation
    std::cout << "\nRunning cellular simulation..." << std::endl;
    cell_system->runTimeCourse(params.total_simulation_time, params.output_interval);
    
    // Analyze results
    auto time_course = cell_system->getTimeCourseData();
    auto metabolic_fluxes = cell_system->getMetabolicFluxes();
    auto signaling_activities = cell_system->getSignalingActivities();
    
    std::cout << "\n=== Simulation Results ===" << std::endl;
    
    // Display some key metabolites over time
    std::cout << "\nKey Metabolite Concentrations (final values):" << std::endl;
    std::cout << std::fixed << std::setprecision(3);
    
    if (time_course.find("glucose") != time_course.end()) {
        auto& glucose_data = time_course.at("glucose");
        std::cout << "Glucose: " << glucose_data.back() << " mM" << std::endl;
    }
    
    if (time_course.find("ATP") != time_course.end()) {
        auto& atp_data = time_course.at("ATP");
        std::cout << "ATP: " << atp_data.back() << " mM" << std::endl;
    }
    
    if (time_course.find("ADP") != time_course.end()) {
        auto& adp_data = time_course.at("ADP");
        std::cout << "ADP: " << adp_data.back() << " mM" << std::endl;
    }
    
    // Display metabolic fluxes
    std::cout << "\nMetabolic Fluxes:" << std::endl;
    for (const auto& [reaction, flux] : metabolic_fluxes) {
        std::cout << reaction << ": " << flux << " mM/s" << std::endl;
    }
    
    // Display signaling activities
    std::cout << "\nSignaling Pathway Activities:" << std::endl;
    for (const auto& [pathway, activity] : signaling_activities) {
        std::cout << pathway << ": " << activity << std::endl;
    }
    
    // Export results
    std::cout << "\nExporting results..." << std::endl;
    cell_system->exportResults("cellular_simulation_results.csv", "csv");
    cell_system->generateReport("cellular_simulation_report.txt");
    
    std::cout << "Results exported to files." << std::endl;
}

/**
 * Demonstrate hepatocyte-specific simulation
 */
void demonstrate_hepatocyte_simulation() {
    std::cout << "\n=== Hepatocyte Simulation Demonstration ===" << std::endl;
    
    auto hepatocyte = CellularSimulationFactory::createHumanHepatocyte();
    
    // Set up metabolic load simulation
    std::map<std::string, double> fed_state = {
        {"glucose", 10.0},     // High glucose (fed state)
        {"fatty_acids", 0.5},  // Low fatty acids
        {"amino_acids", 2.0},  // Moderate amino acids
        {"insulin", 1e-8}      // High insulin (10 nM)
    };
    
    hepatocyte->setExtracellularConditions(fed_state);
    
    // Simulate fed-to-fasted transition
    std::cout << "Simulating fed-to-fasted metabolic transition..." << std::endl;
    
    // Start in fed state, then switch to fasted
    hepatocyte->addMetabolicPerturbation("glucose", -7.0);  // Drop glucose at t=30s
    hepatocyte->addMetabolicPerturbation("fatty_acids", 1.5); // Increase fatty acids
    
    SimulationParameters params;
    params.total_simulation_time = 120.0;  // 2 minutes
    params.time_step = 0.1;
    params.output_interval = 2.0;
    
    hepatocyte->setSimulationParameters(params);
    hepatocyte->runTimeCourse(params.total_simulation_time, params.output_interval);
    
    auto fluxes = hepatocyte->getMetabolicFluxes();
    std::cout << "\nHepatocyte metabolic response:" << std::endl;
    std::cout << "Glycolysis flux: " << (fluxes.count("glycolysis") ? fluxes.at("glycolysis") : 0.0) << std::endl;
    std::cout << "Gluconeogenesis flux: " << (fluxes.count("gluconeogenesis") ? fluxes.at("gluconeogenesis") : 0.0) << std::endl;
    std::cout << "Fatty acid oxidation: " << (fluxes.count("beta_oxidation") ? fluxes.at("beta_oxidation") : 0.0) << std::endl;
}

/**
 * Demonstrate mathematical framework export
 */
void demonstrate_mathematical_framework() {
    std::cout << "\n=== Mathematical Framework Export ===" << std::endl;
    
    auto cell_system = CellularSimulationFactory::createGenericHumanCell();
    
    // Export the mathematical frameworks
    std::cout << "Exporting mathematical frameworks..." << std::endl;
    
    cell_system->exportMathematicalFramework("cellular_mathematical_framework.txt");
    cell_system->exportPSystemRules("p_system_rules.txt");
    cell_system->exportAtomSpaceStructure("atomspace_structure.txt");
    
    std::cout << "Mathematical framework files exported:" << std::endl;
    std::cout << "- cellular_mathematical_framework.txt" << std::endl;
    std::cout << "- p_system_rules.txt" << std::endl;
    std::cout << "- atomspace_structure.txt" << std::endl;
    
    // Display some basic model statistics
    auto molecules = cell_system->listAllMolecules();
    auto reactions = cell_system->listAllReactions();
    auto pathways = cell_system->listAllPathways();
    auto compartments = cell_system->getCompartmentVolumes();
    
    std::cout << "\nModel Statistics:" << std::endl;
    std::cout << "Total molecules: " << molecules.size() << std::endl;
    std::cout << "Total reactions: " << reactions.size() << std::endl;
    std::cout << "Total pathways: " << pathways.size() << std::endl;
    std::cout << "Total compartments: " << compartments.size() << std::endl;
    
    std::cout << "\nCompartment volumes:" << std::endl;
    for (const auto& [compartment, volume] : compartments) {
        std::cout << "  " << compartment << ": " << volume << " pL" << std::endl;
    }
}

/**
 * Demonstrate validation against experimental data
 */
void demonstrate_experimental_validation() {
    std::cout << "\n=== Experimental Validation Demonstration ===" << std::endl;
    
    auto cell_system = CellularSimulationFactory::createGenericHumanCell();
    
    // Create mock experimental data for glucose uptake
    ExperimentalData glucose_uptake_data;
    glucose_uptake_data.experiment_name = "glucose_uptake_kinetics";
    glucose_uptake_data.cell_type = "HEK293";
    glucose_uptake_data.measured_variables = {"glucose", "ATP"};
    
    // Mock time points (0 to 60 seconds)
    for (int i = 0; i <= 60; i += 5) {
        glucose_uptake_data.time_points.push_back(i);
    }
    
    // Mock experimental data with some realistic values
    std::vector<double> glucose_exp = {5.0, 4.8, 4.5, 4.2, 3.9, 3.7, 3.5, 3.3, 3.1, 3.0, 2.9, 2.8, 2.7};
    std::vector<double> atp_exp = {2.0, 2.1, 2.3, 2.4, 2.5, 2.6, 2.7, 2.7, 2.8, 2.8, 2.8, 2.8, 2.8};
    
    glucose_uptake_data.time_series_data.push_back(glucose_exp);
    glucose_uptake_data.time_series_data.push_back(atp_exp);
    
    glucose_uptake_data.experimental_conditions = {
        {"temperature", 310.15},
        {"pH", 7.4},
        {"initial_glucose", 5.0}
    };
    
    // Add experimental data to system
    cell_system->addExperimentalData(glucose_uptake_data);
    
    // Run simulation with same conditions
    cell_system->setExtracellularConditions(glucose_uptake_data.experimental_conditions);
    
    SimulationParameters params;
    params.total_simulation_time = 60.0;
    params.output_interval = 5.0;
    cell_system->setSimulationParameters(params);
    cell_system->runTimeCourse(params.total_simulation_time, params.output_interval);
    
    // Validate against experimental data
    double validation_score = cell_system->validateAgainstExperiment("glucose_uptake_kinetics");
    auto fit_metrics = cell_system->calculateModelFitMetrics();
    
    std::cout << "Validation against experimental glucose uptake data:" << std::endl;
    std::cout << "Overall validation score: " << validation_score << std::endl;
    
    for (const auto& [metric, value] : fit_metrics) {
        std::cout << metric << ": " << value << std::endl;
    }
}

int main() {
    std::cout << "Comprehensive Cellular Simulation System" << std::endl;
    std::cout << "Using OpenCog AtomSpace and Membrane Computing (P-systems)" << std::endl;
    std::cout << "=========================================================" << std::endl;
    
    try {
        // Run demonstrations
        demonstrate_cellular_simulation();
        demonstrate_hepatocyte_simulation();
        demonstrate_mathematical_framework();
        demonstrate_experimental_validation();
        
        std::cout << "\n=== Simulation Complete ===" << std::endl;
        std::cout << "All demonstrations completed successfully!" << std::endl;
        std::cout << "\nThe cellular simulation system includes:" << std::endl;
        std::cout << "✓ OpenCog-inspired AtomSpace for molecular representations" << std::endl;
        std::cout << "✓ Hierarchical P-system rules for membrane interactions" << std::endl;
        std::cout << "✓ Detailed cellular component modeling (membranes, organelles)" << std::endl;
        std::cout << "✓ Comprehensive signaling pathway simulations" << std::endl;
        std::cout << "✓ Specialized mitochondrial implementation" << std::endl;
        std::cout << "✓ Temporal dynamics and reaction kinetics" << std::endl;
        std::cout << "✓ Membrane potential calculations" << std::endl;
        std::cout << "✓ Protein-protein interaction modeling" << std::endl;
        std::cout << "✓ Compartmentalization effects" << std::endl;
        std::cout << "✓ Experimental validation framework" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error during simulation: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}