#ifndef CELLULAR_KINETICS_ENGINE_H
#define CELLULAR_KINETICS_ENGINE_H

#include "../core/atomspace.h"
#include "../core/psystem.h"
#include <vector>
#include <functional>
#include <random>

namespace cellular {

/**
 * Kinetic rate law types
 */
enum class RateLawType {
    MASS_ACTION,           // v = k * [A]^a * [B]^b
    MICHAELIS_MENTEN,      // v = Vmax * [S] / (Km + [S])
    COMPETITIVE_INHIBITION, // v = Vmax * [S] / (Km * (1 + [I]/Ki) + [S])
    ALLOSTERIC,            // Hill equation with cooperativity
    ENZYME_INHIBITION,     // Various inhibition types
    CUSTOM                 // User-defined rate function
};

/**
 * Reaction mechanism types
 */
enum class ReactionMechanism {
    ELEMENTARY,            // Single step reaction
    SEQUENTIAL,            // Ordered sequential mechanism
    RANDOM,                // Random order mechanism
    PING_PONG,             // Ping-pong mechanism
    COOPERATIVE            // Cooperative binding
};

/**
 * Kinetic parameter structure
 */
struct KineticParameters {
    double kcat = 1.0;         // Catalytic rate constant (s⁻¹)
    double km = 1.0;           // Michaelis constant (M)
    double ki = 1.0;           // Inhibition constant (M)
    double kd = 1.0;           // Dissociation constant (M)
    double hill_coefficient = 1.0; // Hill coefficient
    double vmax = 1.0;         // Maximum velocity (M/s)
    double kf = 1.0;           // Forward rate constant
    double kr = 1.0;           // Reverse rate constant
    
    // Temperature dependence (Arrhenius)
    double activation_energy = 50000.0; // J/mol
    double temperature = 310.15;        // K (37°C)
    
    // pH dependence
    double pk1 = 7.0;          // First ionization constant
    double pk2 = 9.0;          // Second ionization constant
    double optimal_ph = 7.4;   // Optimal pH
};

/**
 * Reaction kinetics model
 */
class ReactionKinetics {
public:
    ReactionKinetics(const std::string& reaction_name, RateLawType rate_law);
    
    const std::string& getName() const { return reaction_name_; }
    RateLawType getRateLawType() const { return rate_law_type_; }
    
    // Parameter setting
    void setKineticParameters(const KineticParameters& params);
    const KineticParameters& getKineticParameters() const { return parameters_; }
    
    // Reaction components
    void addReactant(const std::string& reactant, double stoichiometry = 1.0);
    void addProduct(const std::string& product, double stoichiometry = 1.0);
    void addModifier(const std::string& modifier, double effect_strength = 1.0);
    
    // Rate calculation
    double calculateRate(const std::map<std::string, double>& concentrations,
                        double temperature = 310.15, double ph = 7.4) const;
    
    // Custom rate function
    void setCustomRateFunction(std::function<double(const std::map<std::string, double>&)> func);
    
    // Equilibrium
    double calculateEquilibriumConstant(double temperature = 310.15) const;
    bool isAtEquilibrium(const std::map<std::string, double>& concentrations, 
                        double tolerance = 1e-6) const;
    
private:
    std::string reaction_name_;
    RateLawType rate_law_type_;
    KineticParameters parameters_;
    
    std::map<std::string, double> reactants_;
    std::map<std::string, double> products_;
    std::map<std::string, double> modifiers_;
    
    std::function<double(const std::map<std::string, double>&)> custom_rate_function_;
    
    // Rate law implementations
    double massActionRate(const std::map<std::string, double>& concentrations) const;
    double michaelisMentenRate(const std::map<std::string, double>& concentrations) const;
    double competitiveInhibitionRate(const std::map<std::string, double>& concentrations) const;
    double allostericRate(const std::map<std::string, double>& concentrations) const;
    
    // Temperature and pH effects
    double temperatureEffect(double temperature) const;
    double pHEffect(double ph) const;
};

/**
 * Temporal dynamics simulator
 */
class TemporalDynamics {
public:
    enum class IntegrationMethod {
        EULER,              // Simple Euler method
        RUNGE_KUTTA_4,      // Fourth-order Runge-Kutta
        ADAPTIVE_STEP,      // Adaptive step size
        STOCHASTIC,         // Gillespie algorithm
        HYBRID             // Hybrid deterministic/stochastic
    };
    
    TemporalDynamics(IntegrationMethod method = IntegrationMethod::RUNGE_KUTTA_4);
    
    // Integration parameters
    void setIntegrationMethod(IntegrationMethod method) { integration_method_ = method; }
    void setTimeStep(double dt) { time_step_ = dt; }
    void setTolerance(double tolerance) { tolerance_ = tolerance; }
    
    // System state
    void setState(const std::map<std::string, double>& initial_concentrations);
    const std::map<std::string, double>& getState() const { return current_state_; }
    
    // Reaction system
    void addReaction(const ReactionKinetics& reaction);
    void removeReaction(const std::string& reaction_name);
    
    // Simulation
    void step(double dt);
    void simulate(double total_time, double output_interval = -1);
    
    // Stochastic simulation (Gillespie algorithm)
    void gillespieStep();
    void gillespieSimulation(double total_time, int max_steps = 1000000);
    
    // Adaptive time stepping
    void adaptiveStep(double& dt, double tolerance);
    
    // Sensitivity analysis
    std::map<std::string, double> calculateSensitivity(const std::string& parameter,
                                                       double perturbation = 0.01);
    
    // Steady state analysis
    bool findSteadyState(double tolerance = 1e-6, int max_iterations = 10000);
    bool isAtSteadyState(double tolerance = 1e-6) const;
    
    // Stability analysis
    std::vector<std::complex<double>> calculateEigenvalues() const;
    bool isStable() const;
    
    // Time series data
    const std::vector<std::map<std::string, double>>& getTimeSeriesData() const { return time_series_; }
    const std::vector<double>& getTimePoints() const { return time_points_; }
    
    void clearTimeSeriesData();
    
private:
    IntegrationMethod integration_method_;
    double time_step_ = 0.01;
    double tolerance_ = 1e-6;
    double current_time_ = 0.0;
    
    std::map<std::string, double> current_state_;
    std::vector<ReactionKinetics> reactions_;
    
    // Time series storage
    std::vector<std::map<std::string, double>> time_series_;
    std::vector<double> time_points_;
    
    // Random number generation for stochastic simulation
    std::mt19937 rng_;
    std::uniform_real_distribution<double> uniform_dist_;
    
    // Integration methods
    void eulerStep(double dt);
    void rungeKutta4Step(double dt);
    std::map<std::string, double> calculateDerivatives(const std::map<std::string, double>& state) const;
    
    // Stochastic helpers
    double calculateTotalPropensity(const std::map<std::string, double>& state) const;
    int selectReaction(const std::vector<double>& propensities, double random_value) const;
    
    // Jacobian matrix for stability analysis
    std::vector<std::vector<double>> calculateJacobian() const;
};

/**
 * Multi-scale temporal modeling
 */
class MultiScaleTemporalModel {
public:
    MultiScaleTemporalModel();
    
    // Time scale separation
    void addFastTimeScale(double characteristic_time, const std::vector<std::string>& fast_variables);
    void addSlowTimeScale(double characteristic_time, const std::vector<std::string>& slow_variables);
    
    // Fast-slow decomposition
    void performQuasiSteadyStateApproximation(const std::vector<std::string>& fast_variables);
    void performAdiabaticeElimination(const std::vector<std::string>& fast_variables);
    
    // Multi-scale simulation
    void simulateMultiScale(double total_time, double fast_dt, double slow_dt);
    
    // Scale analysis
    std::map<std::string, double> calculateTimeScales() const;
    double getStiffnessRatio() const;
    
private:
    struct TimeScale {
        double characteristic_time;
        std::vector<std::string> variables;
        std::unique_ptr<TemporalDynamics> dynamics;
    };
    
    std::vector<TimeScale> time_scales_;
    double stiffness_ratio_ = 1.0;
};

/**
 * Fluctuation and noise modeling
 */
class FluctuationModel {
public:
    enum class NoiseType {
        WHITE_NOISE,        // Uncorrelated noise
        COLORED_NOISE,      // Correlated noise
        INTRINSIC_NOISE,    // Molecular fluctuations
        EXTRINSIC_NOISE     // Environmental fluctuations
    };
    
    FluctuationModel();
    
    // Noise sources
    void addNoiseSource(const std::string& variable, NoiseType type, double strength);
    void setCorrelationTime(double tau) { correlation_time_ = tau; }
    
    // Fluctuation-dissipation theorem
    void setDiffusionCoefficient(const std::string& variable, double diffusion_coeff);
    
    // Stochastic differential equations
    void addSDEterm(const std::string& variable, 
                   std::function<double(double, double)> drift,
                   std::function<double(double, double)> diffusion);
    
    // Simulation
    void simulateWithFluctuations(TemporalDynamics& dynamics, double dt);
    
    // Analysis
    double calculateVariance(const std::string& variable) const;
    double calculateCorrelationFunction(const std::string& variable, double lag_time) const;
    
private:
    struct NoiseSource {
        NoiseType type;
        double strength;
        double current_value;
    };
    
    std::map<std::string, NoiseSource> noise_sources_;
    std::map<std::string, double> diffusion_coefficients_;
    double correlation_time_ = 1.0;
    
    std::mt19937 rng_;
    std::normal_distribution<double> normal_dist_;
    
    double generateNoise(NoiseType type, double strength, double dt);
    double generateColoredNoise(double current_value, double strength, double dt);
};

/**
 * Comprehensive kinetics engine
 */
class KineticsEngine {
public:
    KineticsEngine();
    
    // Temporal dynamics
    TemporalDynamics& getTemporalDynamics() { return temporal_dynamics_; }
    const TemporalDynamics& getTemporalDynamics() const { return temporal_dynamics_; }
    
    // Multi-scale modeling
    MultiScaleTemporalModel& getMultiScaleModel() { return multiscale_model_; }
    
    // Fluctuations
    FluctuationModel& getFluctuationModel() { return fluctuation_model_; }
    
    // Integration with P-system
    void syncWithPSystem(PSystem& psystem, const AtomSpace& atomspace);
    void updatePSystemFromKinetics(PSystem& psystem, AtomSpace& atomspace);
    
    // Kinetic analysis
    void performKineticAnalysis(const std::string& output_file = "");
    std::map<std::string, double> calculateFluxes() const;
    std::map<std::string, double> calculateElasticities() const;
    
    // Model validation
    void validateAgainstExperimentalData(const std::vector<std::pair<double, std::map<std::string, double>>>& data);
    double calculateRSquared(const std::string& variable, 
                           const std::vector<std::pair<double, double>>& experimental_data) const;
    
    // Parameter estimation
    void fitParameters(const std::vector<std::pair<double, std::map<std::string, double>>>& data,
                      const std::vector<std::string>& parameters_to_fit);
    
    // Simulation control
    void runFullSimulation(double total_time, AtomSpace& atomspace, PSystem& psystem);
    
private:
    TemporalDynamics temporal_dynamics_;
    MultiScaleTemporalModel multiscale_model_;
    FluctuationModel fluctuation_model_;
    
    // Analysis results
    std::map<std::string, double> flux_analysis_;
    std::map<std::string, double> elasticity_analysis_;
    std::map<std::string, double> fitted_parameters_;
    
    void updateAnalysis();
    double objectiveFunction(const std::vector<double>& parameters,
                           const std::vector<std::pair<double, std::map<std::string, double>>>& data);
};

} // namespace cellular

#endif // CELLULAR_KINETICS_ENGINE_H