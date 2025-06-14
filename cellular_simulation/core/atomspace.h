#ifndef CELLULAR_ATOMSPACE_H
#define CELLULAR_ATOMSPACE_H

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

namespace cellular {

// Forward declarations
class Atom;
class AtomSpace;

using AtomPtr = std::shared_ptr<Atom>;
using AtomPtrVector = std::vector<AtomPtr>;

/**
 * Atom types for cellular components
 */
enum class AtomType {
    CONCEPT,        // Generic concept
    MOLECULE,       // Molecular entity
    PROTEIN,        // Protein
    ENZYME,         // Enzyme
    RECEPTOR,       // Cell receptor
    CHANNEL,        // Ion channel
    TRANSPORTER,    // Membrane transporter
    ORGANELLE,      // Cellular organelle
    MEMBRANE,       // Membrane structure
    REACTION,       // Chemical reaction
    PATHWAY,        // Signaling pathway
    COMPLEX,        // Protein complex
    LINK            // Relationship link
};

/**
 * Atom truth value for probabilistic reasoning
 */
struct TruthValue {
    double strength = 1.0;  // [0,1] - confidence in the atom
    double confidence = 1.0; // [0,1] - confidence in the strength
    
    TruthValue() = default;
    TruthValue(double s, double c) : strength(s), confidence(c) {}
};

/**
 * Basic atom in the cellular AtomSpace
 */
class Atom {
public:
    Atom(AtomType type, const std::string& name);
    Atom(AtomType type, const std::string& name, const AtomPtrVector& outgoing);
    
    AtomType getType() const { return type_; }
    const std::string& getName() const { return name_; }
    const AtomPtrVector& getOutgoing() const { return outgoing_; }
    
    const TruthValue& getTruthValue() const { return truth_value_; }
    void setTruthValue(const TruthValue& tv) { truth_value_ = tv; }
    
    // Concentration for molecular entities (mol/L)
    double getConcentration() const { return concentration_; }
    void setConcentration(double conc) { concentration_ = conc; }
    
    // Activity level for reactions/pathways [0,1]
    double getActivity() const { return activity_; }
    void setActivity(double act) { activity_ = act; }
    
    std::string toString() const;
    
private:
    AtomType type_;
    std::string name_;
    AtomPtrVector outgoing_;
    TruthValue truth_value_;
    double concentration_ = 0.0;
    double activity_ = 0.0;
};

/**
 * OpenCog-inspired AtomSpace for cellular components
 */
class AtomSpace {
public:
    AtomSpace() = default;
    
    // Atom creation and management
    AtomPtr addAtom(AtomType type, const std::string& name);
    AtomPtr addAtom(AtomType type, const std::string& name, const AtomPtrVector& outgoing);
    
    AtomPtr findAtom(const std::string& name) const;
    AtomPtrVector findAtomsOfType(AtomType type) const;
    
    // Pattern matching
    AtomPtrVector getIncomingLinks(AtomPtr atom) const;
    AtomPtrVector getOutgoingLinks(AtomPtr atom) const;
    
    // Molecular operations
    void setMolecularConcentration(const std::string& molecule, double concentration);
    double getMolecularConcentration(const std::string& molecule) const;
    
    // Statistics
    size_t getAtomCount() const { return atoms_.size(); }
    size_t getAtomCountByType(AtomType type) const;
    
    void clear();
    
private:
    std::vector<AtomPtr> atoms_;
    std::unordered_map<std::string, AtomPtr> name_index_;
    std::unordered_map<AtomType, std::unordered_set<AtomPtr>> type_index_;
    std::unordered_map<AtomPtr, std::unordered_set<AtomPtr>> incoming_index_;
};

} // namespace cellular

#endif // CELLULAR_ATOMSPACE_H