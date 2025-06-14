#include "psystem.h"
#include <cmath>
#include <algorithm>

namespace cellular {

// ReactionRule implementation
ReactionRule::ReactionRule(const std::string& n, 
                          const std::vector<std::string>& reacts,
                          const std::vector<std::string>& prods,
                          double rate,
                          MembraneType mem) 
    : name(n), reactants(reacts), products(prods), rate_constant(rate), membrane(mem) {
    
    // Initialize stoichiometric coefficients to 1.0
    reactant_stoich.resize(reactants.size(), 1.0);
    product_stoich.resize(products.size(), 1.0);
}

// TransportRule implementation
TransportRule::TransportRule(const std::string& n,
                            const std::string& mol,
                            MembraneType src,
                            MembraneType tgt,
                            double rate)
    : name(n), molecule(mol), source_membrane(src), target_membrane(tgt), transport_rate(rate) {
}

// Membrane implementation
Membrane::Membrane(MembraneType type, const std::string& name) 
    : type_(type), name_(name) {
}

void Membrane::setConcentration(const std::string& molecule, double conc) {
    concentrations_[molecule] = std::max(0.0, conc);
}

double Membrane::getConcentration(const std::string& molecule) const {
    auto it = concentrations_.find(molecule);
    return (it != concentrations_.end()) ? it->second : 0.0;
}

void Membrane::setPermeability(const std::string& molecule, double perm) {
    permeabilities_[molecule] = std::max(0.0, perm);
}

double Membrane::getPermeability(const std::string& molecule) const {
    auto it = permeabilities_.find(molecule);
    return (it != permeabilities_.end()) ? it->second : 0.0;
}

void Membrane::addReactionRule(const ReactionRule& rule) {
    reaction_rules_.push_back(rule);
}

void Membrane::addTransportRule(const TransportRule& rule) {
    transport_rules_.push_back(rule);
}

// PSystem implementation
PSystem::PSystem() {
    // Create standard cellular membranes
    addMembrane(MembraneType::CELL_MEMBRANE, "cell_membrane");
    addMembrane(MembraneType::CYTOPLASM, "cytoplasm");
    addMembrane(MembraneType::NUCLEUS, "nucleus");
    addMembrane(MembraneType::NUCLEAR_ENVELOPE, "nuclear_envelope");
    addMembrane(MembraneType::ER_MEMBRANE, "er_membrane");
    addMembrane(MembraneType::ER_LUMEN, "er_lumen");
    addMembrane(MembraneType::MITOCHONDRIAL_OUTER, "mito_outer");
    addMembrane(MembraneType::MITOCHONDRIAL_INNER, "mito_inner");
    addMembrane(MembraneType::MITOCHONDRIAL_MATRIX, "mito_matrix");
    addMembrane(MembraneType::INTERMEMBRANE_SPACE, "intermembrane_space");
    addMembrane(MembraneType::GOLGI_MEMBRANE, "golgi");
    addMembrane(MembraneType::LYSOSOMAL_MEMBRANE, "lysosome");
    addMembrane(MembraneType::PEROXISOMAL_MEMBRANE, "peroxisome");
    
    // Set typical volumes (in picoliters)
    getMembrane(MembraneType::CYTOPLASM)->setVolume(1000.0);
    getMembrane(MembraneType::NUCLEUS)->setVolume(300.0);
    getMembrane(MembraneType::ER_LUMEN)->setVolume(100.0);
    getMembrane(MembraneType::MITOCHONDRIAL_MATRIX)->setVolume(50.0);
    getMembrane(MembraneType::INTERMEMBRANE_SPACE)->setVolume(10.0);
}

void PSystem::addMembrane(MembraneType type, const std::string& name) {
    membranes_[type] = std::make_unique<Membrane>(type, name);
}

Membrane* PSystem::getMembrane(MembraneType type) {
    auto it = membranes_.find(type);
    return (it != membranes_.end()) ? it->second.get() : nullptr;
}

const Membrane* PSystem::getMembrane(MembraneType type) const {
    auto it = membranes_.find(type);
    return (it != membranes_.end()) ? it->second.get() : nullptr;
}

void PSystem::addReactionRule(const ReactionRule& rule) {
    auto* membrane = getMembrane(rule.membrane);
    if (membrane) {
        membrane->addReactionRule(rule);
    }
}

void PSystem::addTransportRule(const TransportRule& rule) {
    auto* source = getMembrane(rule.source_membrane);
    auto* target = getMembrane(rule.target_membrane);
    if (source && target) {
        source->addTransportRule(rule);
    }
}

void PSystem::step(double dt, AtomSpace& atomspace) {
    // Execute reactions in all membranes
    executeReactions(dt, atomspace);
    
    // Execute transport between membranes
    executeTransport(dt, atomspace);
    
    // Update membrane potentials
    updateMembranePotentials(atomspace);
    
    // Sync concentrations back to atomspace
    updateAtomSpaceConcentrations(atomspace);
}

void PSystem::simulate(double total_time, double dt, AtomSpace& atomspace) {
    double current_time = 0.0;
    while (current_time < total_time) {
        step(dt, atomspace);
        current_time += dt;
    }
}

void PSystem::setInitialConcentration(MembraneType membrane, 
                                     const std::string& molecule, 
                                     double concentration) {
    auto* mem = getMembrane(membrane);
    if (mem) {
        mem->setConcentration(molecule, concentration);
    }
}

void PSystem::executeReactions(double dt, AtomSpace& atomspace) {
    for (auto& [type, membrane] : membranes_) {
        for (const auto& rule : membrane->getReactionRules()) {
            double rate = calculateReactionRate(rule, atomspace);
            
            // Calculate reaction extent
            double extent = rate * dt;
            
            // Update reactant concentrations
            for (size_t i = 0; i < rule.reactants.size(); ++i) {
                double current = membrane->getConcentration(rule.reactants[i]);
                double new_conc = current - extent * rule.reactant_stoich[i];
                membrane->setConcentration(rule.reactants[i], std::max(0.0, new_conc));
            }
            
            // Update product concentrations
            for (size_t i = 0; i < rule.products.size(); ++i) {
                double current = membrane->getConcentration(rule.products[i]);
                double new_conc = current + extent * rule.product_stoich[i];
                membrane->setConcentration(rule.products[i], new_conc);
            }
        }
    }
}

void PSystem::executeTransport(double dt, AtomSpace& atomspace) {
    for (auto& [type, membrane] : membranes_) {
        for (const auto& rule : membrane->getTransportRules()) {
            double rate = calculateTransportRate(rule, atomspace);
            
            auto* source = getMembrane(rule.source_membrane);
            auto* target = getMembrane(rule.target_membrane);
            
            if (source && target) {
                double source_conc = source->getConcentration(rule.molecule);
                double transport_amount = rate * source_conc * dt;
                
                // Ensure we don't transport more than available
                transport_amount = std::min(transport_amount, source_conc);
                
                source->setConcentration(rule.molecule, source_conc - transport_amount);
                target->setConcentration(rule.molecule, 
                    target->getConcentration(rule.molecule) + transport_amount);
            }
        }
    }
}

double PSystem::calculateReactionRate(const ReactionRule& rule, const AtomSpace& atomspace) {
    if (rule.rate_function) {
        return rule.rate_function(atomspace);
    }
    
    // Simple mass action kinetics
    double rate = rule.rate_constant;
    auto* membrane = getMembrane(rule.membrane);
    
    for (size_t i = 0; i < rule.reactants.size(); ++i) {
        double conc = membrane->getConcentration(rule.reactants[i]);
        rate *= std::pow(conc, rule.reactant_stoich[i]);
    }
    
    return rate;
}

double PSystem::calculateTransportRate(const TransportRule& rule, const AtomSpace& atomspace) {
    if (rule.transport_function) {
        return rule.transport_function(atomspace);
    }
    
    return rule.transport_rate;
}

void PSystem::updateAtomSpaceConcentrations(AtomSpace& atomspace) {
    for (auto& [type, membrane] : membranes_) {
        // This is a simplified sync - in reality we'd need more sophisticated handling
        // of compartmentalization in the atomspace
    }
}

void PSystem::updateMembranePotentials(const AtomSpace& atomspace) {
    // Simplified membrane potential calculation
    // In reality, this would use Nernst equation and Goldman-Hodgkin-Katz equation
    
    auto* cytoplasm = getMembrane(MembraneType::CYTOPLASM);
    auto* cell_membrane = getMembrane(MembraneType::CELL_MEMBRANE);
    
    if (cytoplasm && cell_membrane) {
        // Typical resting potential calculation (simplified)
        double k_in = cytoplasm->getConcentration("K+");
        double k_out = 5.0; // Typical extracellular K+ (mM)
        double na_in = cytoplasm->getConcentration("Na+");
        double na_out = 145.0; // Typical extracellular Na+ (mM)
        
        if (k_in > 0 && na_in > 0) {
            // Goldman-Hodgkin-Katz equation (simplified) - corrected
            double pk_pna = 0.04; // Na+ permeability is much lower than K+
            double potential = -61.54 * std::log((pk_pna * na_out + k_out) / (pk_pna * na_in + k_in));
            cell_membrane->setPotential(potential);
        } else {
            // Default resting potential
            cell_membrane->setPotential(-70.0);
        }
    }
}

double PSystem::getTotalMolecules(const std::string& molecule) const {
    double total = 0.0;
    for (const auto& [type, membrane] : membranes_) {
        total += membrane->getConcentration(molecule) * membrane->getVolume();
    }
    return total;
}

std::map<MembraneType, double> PSystem::getCompartmentDistribution(const std::string& molecule) const {
    std::map<MembraneType, double> distribution;
    for (const auto& [type, membrane] : membranes_) {
        distribution[type] = membrane->getConcentration(molecule);
    }
    return distribution;
}

} // namespace cellular