#include "atomspace.h"
#include <sstream>
#include <algorithm>

namespace cellular {

// Atom implementation
Atom::Atom(AtomType type, const std::string& name) 
    : type_(type), name_(name) {
}

Atom::Atom(AtomType type, const std::string& name, const AtomPtrVector& outgoing)
    : type_(type), name_(name), outgoing_(outgoing) {
}

std::string Atom::toString() const {
    std::ostringstream oss;
    oss << "(" << static_cast<int>(type_) << " \"" << name_ << "\"";
    if (!outgoing_.empty()) {
        oss << " [";
        for (size_t i = 0; i < outgoing_.size(); ++i) {
            if (i > 0) oss << ", ";
            oss << outgoing_[i]->getName();
        }
        oss << "]";
    }
    oss << " tv=" << truth_value_.strength << "," << truth_value_.confidence;
    if (concentration_ > 0) {
        oss << " conc=" << concentration_;
    }
    if (activity_ > 0) {
        oss << " act=" << activity_;
    }
    oss << ")";
    return oss.str();
}

// AtomSpace implementation
AtomPtr AtomSpace::addAtom(AtomType type, const std::string& name) {
    // Check if atom already exists
    auto existing = findAtom(name);
    if (existing) {
        return existing;
    }
    
    auto atom = std::make_shared<Atom>(type, name);
    atoms_.push_back(atom);
    name_index_[name] = atom;
    type_index_[type].insert(atom);
    
    return atom;
}

AtomPtr AtomSpace::addAtom(AtomType type, const std::string& name, const AtomPtrVector& outgoing) {
    // For links, generate unique name if not provided
    std::string link_name = name;
    if (link_name.empty()) {
        std::ostringstream oss;
        oss << "link_" << atoms_.size();
        link_name = oss.str();
    }
    
    auto atom = std::make_shared<Atom>(type, link_name, outgoing);
    atoms_.push_back(atom);
    name_index_[link_name] = atom;
    type_index_[type].insert(atom);
    
    // Update incoming links index
    for (auto& target : outgoing) {
        incoming_index_[target].insert(atom);
    }
    
    return atom;
}

AtomPtr AtomSpace::findAtom(const std::string& name) const {
    auto it = name_index_.find(name);
    return (it != name_index_.end()) ? it->second : nullptr;
}

AtomPtrVector AtomSpace::findAtomsOfType(AtomType type) const {
    AtomPtrVector result;
    auto it = type_index_.find(type);
    if (it != type_index_.end()) {
        for (const auto& atom : it->second) {
            result.push_back(atom);
        }
    }
    return result;
}

AtomPtrVector AtomSpace::getIncomingLinks(AtomPtr atom) const {
    AtomPtrVector result;
    auto it = incoming_index_.find(atom);
    if (it != incoming_index_.end()) {
        for (const auto& link : it->second) {
            result.push_back(link);
        }
    }
    return result;
}

AtomPtrVector AtomSpace::getOutgoingLinks(AtomPtr atom) const {
    return atom->getOutgoing();
}

void AtomSpace::setMolecularConcentration(const std::string& molecule, double concentration) {
    auto atom = findAtom(molecule);
    if (atom) {
        atom->setConcentration(concentration);
    }
}

double AtomSpace::getMolecularConcentration(const std::string& molecule) const {
    auto atom = findAtom(molecule);
    return atom ? atom->getConcentration() : 0.0;
}

size_t AtomSpace::getAtomCountByType(AtomType type) const {
    auto it = type_index_.find(type);
    return (it != type_index_.end()) ? it->second.size() : 0;
}

void AtomSpace::clear() {
    atoms_.clear();
    name_index_.clear();
    type_index_.clear();
    incoming_index_.clear();
}

} // namespace cellular