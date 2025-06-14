#ifndef CELLULAR_PSYSTEM_H
#define CELLULAR_PSYSTEM_H

#include "atomspace.h"
#include <functional>
#include <map>

namespace cellular {

/**
 * P-system membrane types
 */
enum class MembraneType {
    CELL_MEMBRANE,      // Outer cell membrane
    NUCLEAR_ENVELOPE,   // Nuclear membrane
    ER_MEMBRANE,        // Endoplasmic reticulum
    MITOCHONDRIAL_OUTER,// Mitochondrial outer membrane
    MITOCHONDRIAL_INNER,// Mitochondrial inner membrane
    GOLGI_MEMBRANE,     // Golgi apparatus
    LYSOSOMAL_MEMBRANE, // Lysosome membrane
    PEROXISOMAL_MEMBRANE,// Peroxisome membrane
    CYTOPLASM,          // Cytoplasmic region
    NUCLEUS,            // Nuclear region
    ER_LUMEN,           // ER interior
    MITOCHONDRIAL_MATRIX,// Mitochondrial matrix
    INTERMEMBRANE_SPACE // Mitochondrial intermembrane space
};

/**
 * Reaction rule for P-system
 */
struct ReactionRule {
    std::string name;
    std::vector<std::string> reactants;      // Input molecules
    std::vector<std::string> products;       // Output molecules
    std::vector<double> reactant_stoich;     // Stoichiometric coefficients
    std::vector<double> product_stoich;      // Stoichiometric coefficients
    double rate_constant = 1.0;              // Reaction rate
    MembraneType membrane = MembraneType::CYTOPLASM; // Where reaction occurs
    std::function<double(const AtomSpace&)> rate_function; // Custom rate function
    
    ReactionRule() = default;
    ReactionRule(const std::string& n, 
                 const std::vector<std::string>& reacts,
                 const std::vector<std::string>& prods,
                 double rate = 1.0,
                 MembraneType mem = MembraneType::CYTOPLASM);
};

/**
 * Transport rule for membrane transport
 */
struct TransportRule {
    std::string name;
    std::string molecule;                    // Transported molecule
    MembraneType source_membrane;            // Source compartment
    MembraneType target_membrane;            // Target compartment
    double transport_rate = 1.0;             // Transport rate
    std::string transporter;                 // Transporter protein (optional)
    bool energy_dependent = false;           // Requires ATP/energy
    std::function<double(const AtomSpace&)> transport_function; // Custom transport
    
    TransportRule() = default;
    TransportRule(const std::string& n,
                  const std::string& mol,
                  MembraneType src,
                  MembraneType tgt,
                  double rate = 1.0);
};

/**
 * Membrane compartment in P-system
 */
class Membrane {
public:
    Membrane(MembraneType type, const std::string& name);
    
    MembraneType getType() const { return type_; }
    const std::string& getName() const { return name_; }
    
    // Molecular concentrations in this compartment
    void setConcentration(const std::string& molecule, double conc);
    double getConcentration(const std::string& molecule) const;
    
    // Membrane properties
    void setVolume(double volume) { volume_ = volume; }
    double getVolume() const { return volume_; }
    
    void setPotential(double potential) { potential_ = potential; }
    double getPotential() const { return potential_; }
    
    void setPermeability(const std::string& molecule, double perm);
    double getPermeability(const std::string& molecule) const;
    
    // Rules in this membrane
    void addReactionRule(const ReactionRule& rule);
    void addTransportRule(const TransportRule& rule);
    
    const std::vector<ReactionRule>& getReactionRules() const { return reaction_rules_; }
    const std::vector<TransportRule>& getTransportRules() const { return transport_rules_; }
    
private:
    MembraneType type_;
    std::string name_;
    std::map<std::string, double> concentrations_;
    std::map<std::string, double> permeabilities_;
    std::vector<ReactionRule> reaction_rules_;
    std::vector<TransportRule> transport_rules_;
    double volume_ = 1.0;     // Compartment volume (pL)
    double potential_ = 0.0;  // Membrane potential (mV)
};

/**
 * P-system engine for membrane computing
 */
class PSystem {
public:
    PSystem();
    
    // Membrane management
    void addMembrane(MembraneType type, const std::string& name);
    Membrane* getMembrane(MembraneType type);
    const Membrane* getMembrane(MembraneType type) const;
    
    // Rule management
    void addReactionRule(const ReactionRule& rule);
    void addTransportRule(const TransportRule& rule);
    
    // Simulation
    void step(double dt, AtomSpace& atomspace);
    void simulate(double total_time, double dt, AtomSpace& atomspace);
    
    // Initial conditions
    void setInitialConcentration(MembraneType membrane, 
                               const std::string& molecule, 
                               double concentration);
    
    // Membrane potential calculation
    void updateMembranePotentials(const AtomSpace& atomspace);
    
    // Statistics
    double getTotalMolecules(const std::string& molecule) const;
    std::map<MembraneType, double> getCompartmentDistribution(const std::string& molecule) const;
    
private:
    std::map<MembraneType, std::unique_ptr<Membrane>> membranes_;
    
    // Helper functions
    void executeReactions(double dt, AtomSpace& atomspace);
    void executeTransport(double dt, AtomSpace& atomspace);
    double calculateReactionRate(const ReactionRule& rule, const AtomSpace& atomspace);
    double calculateTransportRate(const TransportRule& rule, const AtomSpace& atomspace);
    void updateAtomSpaceConcentrations(AtomSpace& atomspace);
};

} // namespace cellular

#endif // CELLULAR_PSYSTEM_H