#ifndef CELLULAR_MITOCHONDRIA_H
#define CELLULAR_MITOCHONDRIA_H

#include "../core/atomspace.h"
#include "../core/psystem.h"
#include "../organelles/organelles.h"
#include <vector>
#include <map>
#include <complex>

namespace cellular {

/**
 * Mitochondrial DNA model
 */
class MitochondrialDNA {
public:
    MitochondrialDNA();
    
    // mtDNA properties
    void setmtDNACopyNumber(int copies) { mtdna_copies_ = copies; }
    int getmtDNACopyNumber() const { return mtdna_copies_; }
    
    void setGenomeSize(int size) { genome_size_ = size; }
    int getGenomeSize() const { return genome_size_; }
    
    // Gene expression
    void addmtDNAGene(const std::string& gene_name, double expression_rate);
    double getGeneExpression(const std::string& gene_name) const;
    
    // mtDNA replication
    void setReplicationRate(double rate) { replication_rate_ = rate; }
    void simulateReplication(double dt, double cellular_energy_state);
    
    // mtDNA maintenance and repair
    void setMutationRate(double rate) { mutation_rate_ = rate; }
    void addRepairMechanism(const std::string& mechanism, double efficiency);
    void simulateMaintenance(double dt);
    
    // Inheritance patterns
    void simulateMaternelInheritance(double heteroplasmy_level);
    
private:
    int mtdna_copies_ = 1000;           // Typical mtDNA copies per mitochondrion
    int genome_size_ = 16569;           // Human mtDNA size (bp)
    double replication_rate_ = 0.1;     // per hour
    double mutation_rate_ = 1e-5;       // per base per replication
    
    std::map<std::string, double> gene_expression_;
    std::map<std::string, double> repair_mechanisms_;
    double damage_level_ = 0.0;
    double heteroplasmy_ = 0.0;
};

/**
 * Electron Transport Chain Complex
 */
class ETCComplex {
public:
    ETCComplex(const std::string& name, int complex_number);
    
    const std::string& getName() const { return name_; }
    int getComplexNumber() const { return complex_number_; }
    
    // Subunit composition
    void addSubunit(const std::string& subunit_name, bool mtdna_encoded, double activity);
    void setSubunitActivity(const std::string& subunit_name, double activity);
    
    // Electron transport
    void setElectronCarriers(const std::vector<std::string>& carriers);
    void setRedoxPotentials(const std::vector<double>& potentials);
    
    // Proton pumping
    void setProtonStoichiometry(double protons_per_electron);
    double calculateProtonPumping(double electron_flow_rate) const;
    
    // Activity calculation
    double calculateActivity(const std::map<std::string, double>& substrate_concentrations,
                           double membrane_potential, double ph_gradient) const;
    
    // Inhibition and regulation
    void addInhibitor(const std::string& inhibitor, double ki);
    void addRegulator(const std::string& regulator, double effect_strength);
    
private:
    std::string name_;
    int complex_number_;
    std::map<std::string, std::pair<bool, double>> subunits_; // mtDNA_encoded, activity
    std::vector<std::string> electron_carriers_;
    std::vector<double> redox_potentials_;
    double proton_stoichiometry_ = 1.0;
    std::map<std::string, double> inhibitors_;
    std::map<std::string, double> regulators_;
};

/**
 * ATP Synthase (Complex V)
 */
class ATPSynthase {
public:
    ATPSynthase();
    
    // F0 and F1 components
    void setF0Properties(double proton_channel_conductance, int c_subunit_number);
    void setF1Properties(double catalytic_efficiency, double atp_binding_affinity);
    
    // ATP synthesis mechanism
    double calculateATPSynthesis(double proton_gradient, double adp_concentration,
                               double pi_concentration, double membrane_potential) const;
    
    // Rotational mechanism
    void setRotationalProperties(double rotation_rate, double coupling_efficiency);
    double calculateRotationalVelocity(double proton_flow_rate) const;
    
    // Regulation
    void setInhibitorySites(const std::map<std::string, double>& inhibitors);
    void addRegulatedStep(const std::string& regulator, double regulation_strength);
    
private:
    double f0_conductance_ = 1.0;
    int c_subunit_number_ = 12;
    double f1_efficiency_ = 0.9;
    double atp_binding_kd_ = 1e-6;
    double rotation_rate_ = 100.0;  // rotations per second
    double coupling_efficiency_ = 0.95;
    std::map<std::string, double> inhibitors_;
    std::map<std::string, double> regulators_;
};

/**
 * Mitochondrial membrane model
 */
class MitochondrialMembrane {
public:
    enum class MembraneLayer { OUTER, INNER };
    
    MitochondrialMembrane(MembraneLayer layer);
    
    MembraneLayer getLayer() const { return layer_; }
    
    // Membrane properties
    void setPermeability(const std::string& molecule, double permeability);
    double getPermeability(const std::string& molecule) const;
    
    void setSurfaceArea(double area) { surface_area_ = area; }
    double getSurfaceArea() const { return surface_area_; }
    
    // Cristae structure (inner membrane only)
    void setCristaeProperties(double folding_factor, double surface_enhancement);
    double getCristaeEnhancement() const { return cristae_enhancement_; }
    
    // Membrane potential
    void setMembranePotential(double potential) { membrane_potential_ = potential; }
    double getMembranePotential() const { return membrane_potential_; }
    
    // Transport systems
    void addTransporter(const std::string& transporter_name, const std::string& substrate,
                       double vmax, double km, bool electrogenic = false);
    
    void addChannel(const std::string& channel_name, const std::string& ion,
                   double conductance, double selectivity);
    
    // Membrane dynamics
    double calculateTransport(const std::string& substrate, double concentration_gradient,
                            double membrane_potential = 0.0) const;
    
private:
    MembraneLayer layer_;
    double surface_area_ = 1.0;         // μm²
    double cristae_enhancement_ = 3.0;   // Inner membrane folding factor
    double membrane_potential_ = -180.0; // mV (inner membrane)
    
    struct Transporter {
        std::string substrate;
        double vmax;
        double km;
        bool electrogenic;
    };
    
    struct Channel {
        std::string ion;
        double conductance;
        double selectivity;
    };
    
    std::map<std::string, double> permeabilities_;
    std::map<std::string, Transporter> transporters_;
    std::map<std::string, Channel> channels_;
};

/**
 * Comprehensive mitochondrial model
 */
class Mitochondrion : public Organelle {
public:
    Mitochondrion();
    
    // Mitochondrial compartments
    MitochondrialMembrane& getOuterMembrane() { return outer_membrane_; }
    MitochondrialMembrane& getInnerMembrane() { return inner_membrane_; }
    const MitochondrialMembrane& getOuterMembrane() const { return outer_membrane_; }
    const MitochondrialMembrane& getInnerMembrane() const { return inner_membrane_; }
    
    // Compartment volumes
    void setMatrixVolume(double volume) { matrix_volume_ = volume; }
    void setIntermembraneVolume(double volume) { intermembrane_volume_ = volume; }
    double getMatrixVolume() const { return matrix_volume_; }
    double getIntermembraneVolume() const { return intermembrane_volume_; }
    
    // mtDNA operations
    MitochondrialDNA& getmtDNA() { return mtdna_; }
    const MitochondrialDNA& getmtDNA() const { return mtdna_; }
    
    // Electron transport chain
    void addETCComplex(const ETCComplex& complex);
    ETCComplex* getETCComplex(int complex_number);
    
    ATPSynthase& getATPSynthase() { return atp_synthase_; }
    const ATPSynthase& getATPSynthase() const { return atp_synthase_; }
    
    // Respiratory control
    void setRespiratoryState(int state) { respiratory_state_ = state; }
    int getRespiratoryState() const { return respiratory_state_; }
    
    // Substrate oxidation
    void addSubstrateOxidation(const std::string& substrate, const std::string& complex_entry,
                              double oxidation_rate, double nadh_yield, double fadh2_yield);
    
    // TCA cycle
    void initializeTCAcycle();
    void updateTCAcycle(double dt, const std::map<std::string, double>& substrate_concentrations);
    
    // Fatty acid oxidation
    void initializeBetaOxidation();
    void updateBetaOxidation(double dt, const std::map<std::string, double>& fatty_acid_concentrations);
    
    // Calcium handling
    void setCalciumProperties(double uniporter_activity, double efflux_rate, double buffering_capacity);
    double calculateCalciumUptake(double cytosolic_calcium, double membrane_potential) const;
    
    // ROS production and detoxification
    void setROSProperties(double production_rate, double antioxidant_capacity);
    double calculateROSProduction(double electron_transport_rate) const;
    double calculateROSDetoxification(double ros_level) const;
    
    // Mitochondrial dynamics
    void setDynamicsProperties(double fission_rate, double fusion_rate, double motility);
    void updateMitochondrialDynamics(double dt, double cellular_stress_level);
    
    // Bioenergetics calculation
    double calculateMembraneGradient() const;
    double calculateProtonmotiveForce() const;
    double calculateATPProductionRate(const std::map<std::string, double>& substrate_concentrations) const;
    double calculateOxygenConsumption() const;
    
    // Quality control
    void updateMitophagy(double dt, double damage_threshold);
    void updateBiogenesis(double dt, double energy_demand);
    
    void updateOrganelle(double dt, AtomSpace& atomspace, PSystem& psystem) override;
    void initializeOrganelle(AtomSpace& atomspace, PSystem& psystem) override;
    
private:
    MitochondrialMembrane outer_membrane_;
    MitochondrialMembrane inner_membrane_;
    MitochondrialDNA mtdna_;
    ATPSynthase atp_synthase_;
    
    std::vector<ETCComplex> etc_complexes_;
    
    double matrix_volume_ = 0.3;        // pL
    double intermembrane_volume_ = 0.1; // pL
    int respiratory_state_ = 3;         // State 3 (ADP-stimulated)
    
    struct SubstrateOxidation {
        std::string complex_entry;
        double oxidation_rate;
        double nadh_yield;
        double fadh2_yield;
    };
    
    std::map<std::string, SubstrateOxidation> substrate_oxidations_;
    
    // TCA cycle enzymes and rates
    std::map<std::string, double> tca_enzymes_;
    std::map<std::string, double> tca_intermediates_;
    
    // Beta oxidation
    std::map<std::string, double> beta_oxidation_enzymes_;
    
    // Calcium handling
    double calcium_uniporter_activity_ = 1.0;
    double calcium_efflux_rate_ = 0.1;
    double calcium_buffering_capacity_ = 100.0; // μM
    
    // ROS metabolism
    double ros_production_rate_ = 0.02; // % of electron flow
    double antioxidant_capacity_ = 10.0;
    
    // Dynamics
    double fission_rate_ = 0.1;
    double fusion_rate_ = 0.1;
    double motility_rate_ = 1.0;
    
    // Quality control
    double damage_level_ = 0.0;
    double biogenesis_rate_ = 0.01;
    
    void simulateElectronTransport(double dt, const std::map<std::string, double>& reducing_equivalents);
    void simulateATPSynthesis(double dt, double proton_gradient);
    void updateMembraneGradients(double dt);
    void simulateMetabolicFlexibility(double dt, const std::map<std::string, double>& fuel_availability);
};

/**
 * Mitochondrial network model
 */
class MitochondrialNetwork {
public:
    MitochondrialNetwork();
    
    // Network composition
    void addMitochondrion(std::unique_ptr<Mitochondrion> mitochondrion);
    void setNetworkSize(int number_of_mitochondria);
    
    // Network dynamics
    void simulateFissionFusion(double dt, double stress_level);
    void simulateMitochondrialMotility(double dt, double energy_demand_spatial_distribution);
    
    // Collective properties
    double getTotalATPProduction() const;
    double getTotalOxygenConsumption() const;
    double getAverageMembraneGradient() const;
    
    // Quality control at network level
    void simulateNetworkMitophagy(double dt, double damage_threshold);
    void simulateNetworkBiogenesis(double dt, double energy_demand);
    
    // Metabolic flexibility
    void adjustNetworkMetabolism(const std::map<std::string, double>& nutrient_availability);
    
    void updateNetwork(double dt, AtomSpace& atomspace, PSystem& psystem);
    
private:
    std::vector<std::unique_ptr<Mitochondrion>> mitochondria_;
    double network_connectivity_ = 0.7; // Fraction of connected mitochondria
    double spatial_heterogeneity_ = 0.3; // Spatial distribution variability
    
    void calculateNetworkTopology();
    void updateNetworkConnectivity(double dt);
};

} // namespace cellular

#endif // CELLULAR_MITOCHONDRIA_H