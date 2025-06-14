#ifndef CELLULAR_MEMBRANE_DYNAMICS_H
#define CELLULAR_MEMBRANE_DYNAMICS_H

#include "../core/atomspace.h"
#include "../core/psystem.h"
#include <vector>
#include <map>

namespace cellular {

/**
 * Phospholipid types in cell membranes
 */
enum class PhospholipidType {
    PHOSPHATIDYLCHOLINE,     // PC - most abundant
    PHOSPHATIDYLETHANOLAMINE,// PE - second most abundant
    PHOSPHATIDYLSERINE,      // PS - asymmetrically distributed
    PHOSPHATIDYLINOSITOL,    // PI - signaling lipid
    SPHINGOMYELIN,           // SM - outer leaflet
    CHOLESTEROL              // Membrane fluidity regulator
};

/**
 * Membrane transport mechanisms
 */
enum class TransportType {
    PASSIVE_DIFFUSION,       // Simple diffusion
    FACILITATED_DIFFUSION,   // Channel/carrier mediated
    ACTIVE_TRANSPORT,        // ATP-dependent
    SYMPORT,                 // Co-transport same direction
    ANTIPORT,                // Counter-transport
    ENDOCYTOSIS,             // Vesicle formation
    EXOCYTOSIS               // Vesicle fusion
};

/**
 * Ion channel types
 */
enum class ChannelType {
    VOLTAGE_GATED_NA,        // Voltage-gated sodium
    VOLTAGE_GATED_K,         // Voltage-gated potassium
    VOLTAGE_GATED_CA,        // Voltage-gated calcium
    LIGAND_GATED,            // Ligand-gated channels
    MECHANOSENSITIVE,        // Mechanosensitive channels
    LEAK_CHANNEL             // Leak channels
};

/**
 * Membrane receptor types
 */
enum class ReceptorType {
    GPCR,                    // G-protein coupled receptor
    RECEPTOR_TYROSINE_KINASE,// RTK
    ION_CHANNEL_RECEPTOR,    // Ligand-gated ion channel
    ENZYME_LINKED_RECEPTOR,  // Enzyme-linked receptor
    NUCLEAR_RECEPTOR         // Intracellular receptor
};

/**
 * Phospholipid bilayer model
 */
class PhospholipidBilayer {
public:
    PhospholipidBilayer();
    
    // Lipid composition
    void setLipidComposition(PhospholipidType lipid, double outer_leaflet, double inner_leaflet);
    double getLipidFraction(PhospholipidType lipid, bool outer_leaflet = true) const;
    
    // Membrane properties
    double getThickness() const { return thickness_; }
    double getFluidity() const { return fluidity_; }
    double getAsymmetry() const { return asymmetry_; }
    
    void updateProperties();
    
    // Lipid dynamics
    void simulateFlipFlop(double dt);
    void simulateLateralDiffusion(double dt);
    
private:
    std::map<PhospholipidType, double> outer_leaflet_;
    std::map<PhospholipidType, double> inner_leaflet_;
    double thickness_ = 4.0;  // nm
    double fluidity_ = 0.5;   // [0,1]
    double asymmetry_ = 0.1;  // Asymmetry index
};

/**
 * Ion channel model
 */
class IonChannel {
public:
    IonChannel(ChannelType type, const std::string& ion);
    
    ChannelType getType() const { return type_; }
    const std::string& getIon() const { return ion_; }
    
    // Channel properties
    void setMaxConductance(double conductance) { max_conductance_ = conductance; }
    double getMaxConductance() const { return max_conductance_; }
    
    void setOpenProbability(double prob) { open_probability_ = prob; }
    double getOpenProbability() const { return open_probability_; }
    
    // Gating mechanisms
    void updateVoltageGating(double membrane_potential);
    void updateLigandGating(double ligand_concentration);
    
    // Current calculation
    double calculateCurrent(double membrane_potential, double ion_gradient) const;
    
private:
    ChannelType type_;
    std::string ion_;
    double max_conductance_ = 1.0;  // pS (picosiemens)
    double open_probability_ = 0.0; // [0,1]
    double half_activation_ = -50.0; // mV
    double slope_factor_ = 10.0;     // mV
};

/**
 * Membrane transporter model
 */
class MembraneTransporter {
public:
    MembraneTransporter(const std::string& name, TransportType type);
    
    const std::string& getName() const { return name_; }
    TransportType getType() const { return type_; }
    
    // Transport properties
    void setKm(double km) { km_ = km; }           // Michaelis constant
    void setVmax(double vmax) { vmax_ = vmax; }   // Maximum velocity
    
    double getKm() const { return km_; }
    double getVmax() const { return vmax_; }
    
    // Transport rate calculation
    double calculateTransportRate(double substrate_conc, double atp_conc = 0.0) const;
    
    // Stoichiometry
    void setStoichiometry(const std::string& molecule, int coefficient);
    int getStoichiometry(const std::string& molecule) const;
    
private:
    std::string name_;
    TransportType type_;
    double km_ = 1.0;      // mM
    double vmax_ = 1.0;    // molecules/s
    std::map<std::string, int> stoichiometry_;
};

/**
 * Membrane receptor model
 */
class MembraneReceptor {
public:
    MembraneReceptor(const std::string& name, ReceptorType type);
    
    const std::string& getName() const { return name_; }
    ReceptorType getType() const { return type_; }
    
    // Binding properties
    void setKd(double kd) { kd_ = kd; }          // Dissociation constant
    void setMaxBinding(double max) { max_binding_ = max; } // Maximum binding sites
    
    double getKd() const { return kd_; }
    double getMaxBinding() const { return max_binding_; }
    
    // Receptor activation
    double calculateOccupancy(double ligand_concentration) const;
    double calculateActivity(double ligand_concentration) const;
    
    // Downstream signaling
    void setDownstreamTargets(const std::vector<std::string>& targets);
    const std::vector<std::string>& getDownstreamTargets() const { return downstream_targets_; }
    
private:
    std::string name_;
    ReceptorType type_;
    double kd_ = 1e-9;     // M (nM range typical)
    double max_binding_ = 1000.0; // Number of receptors
    std::vector<std::string> downstream_targets_;
    double hill_coefficient_ = 1.0; // Cooperativity
};

/**
 * Comprehensive membrane dynamics system
 */
class MembraneDynamics {
public:
    MembraneDynamics();
    
    // Bilayer management
    void initializeBilayer(MembraneType membrane_type);
    PhospholipidBilayer* getBilayer(MembraneType membrane_type);
    
    // Channel management
    void addIonChannel(MembraneType membrane, const IonChannel& channel);
    std::vector<IonChannel>& getIonChannels(MembraneType membrane);
    
    // Transporter management
    void addTransporter(MembraneType membrane, const MembraneTransporter& transporter);
    std::vector<MembraneTransporter>& getTransporters(MembraneType membrane);
    
    // Receptor management
    void addReceptor(MembraneType membrane, const MembraneReceptor& receptor);
    std::vector<MembraneReceptor>& getReceptors(MembraneType membrane);
    
    // Simulation
    void simulate(double dt, PSystem& psystem, AtomSpace& atomspace);
    
    // Membrane potential calculation
    double calculateMembranePotential(MembraneType membrane, const AtomSpace& atomspace);
    
    // Transport rate calculations
    void updateTransportRates(PSystem& psystem, const AtomSpace& atomspace);
    
    // Receptor signaling
    void updateReceptorSignaling(AtomSpace& atomspace);
    
private:
    std::map<MembraneType, std::unique_ptr<PhospholipidBilayer>> bilayers_;
    std::map<MembraneType, std::vector<IonChannel>> ion_channels_;
    std::map<MembraneType, std::vector<MembraneTransporter>> transporters_;
    std::map<MembraneType, std::vector<MembraneReceptor>> receptors_;
    
    // Helper functions
    void simulateElectrophysiology(double dt, PSystem& psystem, AtomSpace& atomspace);
    void simulateTransport(double dt, PSystem& psystem, AtomSpace& atomspace);
    void simulateReceptorSignaling(double dt, AtomSpace& atomspace);
};

} // namespace cellular

#endif // CELLULAR_MEMBRANE_DYNAMICS_H