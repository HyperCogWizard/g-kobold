#ifndef CELLULAR_ORGANELLES_H
#define CELLULAR_ORGANELLES_H

#include "../core/atomspace.h"
#include "../core/psystem.h"
#include <vector>
#include <map>
#include <memory>

namespace cellular {

/**
 * Base class for all organelles
 */
class Organelle {
public:
    Organelle(const std::string& name, MembraneType membrane_type);
    virtual ~Organelle() = default;
    
    const std::string& getName() const { return name_; }
    MembraneType getMembraneType() const { return membrane_type_; }
    
    // Volume and surface area
    void setVolume(double volume) { volume_ = volume; }
    double getVolume() const { return volume_; }
    
    void setSurfaceArea(double area) { surface_area_ = area; }
    double getSurfaceArea() const { return surface_area_; }
    
    // Molecular concentrations
    void setConcentration(const std::string& molecule, double concentration);
    double getConcentration(const std::string& molecule) const;
    
    // Virtual methods for organelle-specific functions
    virtual void updateOrganelle(double dt, AtomSpace& atomspace, PSystem& psystem) = 0;
    virtual void initializeOrganelle(AtomSpace& atomspace, PSystem& psystem) = 0;
    
protected:
    std::string name_;
    MembraneType membrane_type_;
    double volume_ = 1.0;        // pL
    double surface_area_ = 1.0;  // μm²
    std::map<std::string, double> concentrations_;
};

/**
 * Endoplasmic Reticulum model
 */
class EndoplasmicReticulum : public Organelle {
public:
    EndoplasmicReticulum();
    
    // ER types
    enum class ERType { ROUGH, SMOOTH };
    void setERType(ERType type) { er_type_ = type; }
    ERType getERType() const { return er_type_; }
    
    // Ribosome binding (for rough ER)
    void setRibosomeDensity(double density) { ribosome_density_ = density; }
    double getRibosomeDensity() const { return ribosome_density_; }
    
    // Protein synthesis and folding
    void addProteinSynthesisReaction(const std::string& protein, double synthesis_rate);
    void addProteinFoldingReaction(const std::string& protein, double folding_rate, double misfolding_rate);
    
    // Calcium storage and release
    void setCalciumProperties(double storage_capacity, double release_rate, double uptake_rate);
    double calculateCalciumRelease(double ip3_concentration) const;
    
    // Lipid synthesis (for smooth ER)
    void addLipidSynthesisReaction(const std::string& lipid, double synthesis_rate);
    
    // Drug metabolism (for smooth ER)
    void addCytochromeP450Reaction(const std::string& substrate, double metabolism_rate);
    
    void updateOrganelle(double dt, AtomSpace& atomspace, PSystem& psystem) override;
    void initializeOrganelle(AtomSpace& atomspace, PSystem& psystem) override;
    
private:
    ERType er_type_ = ERType::ROUGH;
    double ribosome_density_ = 100.0; // ribosomes/μm²
    double calcium_storage_capacity_ = 0.5; // mM
    double calcium_release_rate_ = 5.0;
    double calcium_uptake_rate_ = 2.0;
    
    std::map<std::string, double> protein_synthesis_rates_;
    std::map<std::string, std::pair<double, double>> protein_folding_rates_; // folding, misfolding
    std::map<std::string, double> lipid_synthesis_rates_;
    std::map<std::string, double> p450_reactions_;
};

/**
 * Golgi Apparatus model
 */
class GolgiApparatus : public Organelle {
public:
    GolgiApparatus();
    
    // Golgi compartments
    enum class GolgiCompartment { CIS, MEDIAL, TRANS, TGN }; // Trans-Golgi Network
    
    // Protein processing and modification
    void addGlycosylationReaction(const std::string& protein, const std::string& sugar,
                                 GolgiCompartment compartment, double rate);
    
    void addPhosphorylationReaction(const std::string& protein, 
                                   GolgiCompartment compartment, double rate);
    
    void addSulfationReaction(const std::string& protein,
                             GolgiCompartment compartment, double rate);
    
    // Vesicle transport
    void setVesicleTransportRate(GolgiCompartment from, GolgiCompartment to, double rate);
    void addSecretoryVesicleFormation(const std::string& cargo, double formation_rate);
    
    // Lipid modification
    void addLipidModificationReaction(const std::string& lipid, double modification_rate);
    
    void updateOrganelle(double dt, AtomSpace& atomspace, PSystem& psystem) override;
    void initializeOrganelle(AtomSpace& atomspace, PSystem& psystem) override;
    
private:
    struct ProcessingReaction {
        std::string substrate;
        std::string modifier;
        GolgiCompartment compartment;
        double rate;
    };
    
    std::vector<ProcessingReaction> glycosylation_reactions_;
    std::vector<ProcessingReaction> phosphorylation_reactions_;
    std::vector<ProcessingReaction> sulfation_reactions_;
    std::map<std::pair<GolgiCompartment, GolgiCompartment>, double> transport_rates_;
    std::map<std::string, double> secretory_vesicle_rates_;
    std::map<std::string, double> lipid_modification_rates_;
};

/**
 * Lysosome model
 */
class Lysosome : public Organelle {
public:
    Lysosome();
    
    // pH regulation
    void setpH(double ph) { ph_ = ph; }
    double getpH() const { return ph_; }
    
    // Hydrolytic enzymes
    void addHydrolyticEnzyme(const std::string& enzyme, const std::string& substrate,
                            double activity, double ph_optimum);
    
    // Autophagy
    void setAutophagyRate(double rate) { autophagy_rate_ = rate; }
    void addAutophagyTarget(const std::string& target, double degradation_rate);
    
    // Endocytosis and degradation
    void addEndocyticCargo(const std::string& cargo, double uptake_rate, double degradation_rate);
    
    // pH buffering and proton pumping
    void updatepH(double dt, double proton_pump_activity);
    
    void updateOrganelle(double dt, AtomSpace& atomspace, PSystem& psystem) override;
    void initializeOrganelle(AtomSpace& atomspace, PSystem& psystem) override;
    
private:
    double ph_ = 4.5; // Acidic pH
    double autophagy_rate_ = 0.1;
    
    struct HydrolyticEnzyme {
        std::string enzyme_name;
        std::string substrate;
        double activity;
        double ph_optimum;
    };
    
    std::vector<HydrolyticEnzyme> hydrolytic_enzymes_;
    std::map<std::string, double> autophagy_targets_;
    std::map<std::string, std::pair<double, double>> endocytic_cargo_; // uptake, degradation
};

/**
 * Peroxisome model
 */
class Peroxisome : public Organelle {
public:
    Peroxisome();
    
    // Fatty acid oxidation
    void addBetaOxidationReaction(const std::string& fatty_acid, double oxidation_rate);
    
    // ROS metabolism
    void setCatalaseActivity(double activity) { catalase_activity_ = activity; }
    void setPeroxidaseActivity(double activity) { peroxidase_activity_ = activity; }
    
    // Hydrogen peroxide production and consumption
    double calculateH2O2Production(const std::map<std::string, double>& substrate_concentrations) const;
    double calculateH2O2Consumption(double h2o2_concentration) const;
    
    // Glyoxylate cycle (in plants)
    void addGlyoxylateCycleReaction(const std::string& substrate, double rate);
    
    // Plasmalogens synthesis
    void addPlasmalogenSynthesis(const std::string& plasmalogen, double synthesis_rate);
    
    void updateOrganelle(double dt, AtomSpace& atomspace, PSystem& psystem) override;
    void initializeOrganelle(AtomSpace& atomspace, PSystem& psystem) override;
    
private:
    double catalase_activity_ = 1000.0;  // High activity for H2O2 decomposition
    double peroxidase_activity_ = 100.0;
    
    std::map<std::string, double> beta_oxidation_rates_;
    std::map<std::string, double> glyoxylate_cycle_rates_;
    std::map<std::string, double> plasmalogen_synthesis_rates_;
};

/**
 * Nucleus model
 */
class Nucleus : public Organelle {
public:
    Nucleus();
    
    // DNA and chromatin
    void setDNAContent(double content) { dna_content_ = content; }
    double getDNAContent() const { return dna_content_; }
    
    void setChromatinState(const std::string& gene, double accessibility);
    double getChromatinAccessibility(const std::string& gene) const;
    
    // Transcription
    void addTranscriptionFactor(const std::string& tf_name, const std::string& target_gene,
                               double binding_affinity, double activation_strength);
    
    void addGeneTranscription(const std::string& gene, double basal_rate, double enhancement_factor);
    
    // RNA processing
    void addRNAProcessing(const std::string& pre_mrna, const std::string& mature_mrna,
                         double processing_rate);
    
    // Nuclear export
    void addNuclearExport(const std::string& rna_type, double export_rate);
    
    // DNA repair
    void addDNARepairPathway(const std::string& damage_type, double repair_rate);
    
    void updateOrganelle(double dt, AtomSpace& atomspace, PSystem& psystem) override;
    void initializeOrganelle(AtomSpace& atomspace, PSystem& psystem) override;
    
private:
    double dna_content_ = 6.0; // pg (diploid human genome)
    
    struct TranscriptionFactor {
        std::string name;
        std::string target_gene;
        double binding_affinity;
        double activation_strength;
    };
    
    std::map<std::string, double> chromatin_accessibility_;
    std::vector<TranscriptionFactor> transcription_factors_;
    std::map<std::string, std::pair<double, double>> gene_transcription_; // basal_rate, enhancement
    std::map<std::string, std::pair<std::string, double>> rna_processing_; // target, rate
    std::map<std::string, double> nuclear_export_rates_;
    std::map<std::string, double> dna_repair_rates_;
};

/**
 * Cytoskeleton model
 */
class Cytoskeleton : public Organelle {
public:
    Cytoskeleton();
    
    // Cytoskeletal elements
    enum class FilamentType { MICROFILAMENTS, MICROTUBULES, INTERMEDIATE_FILAMENTS };
    
    void setFilamentDensity(FilamentType type, double density);
    double getFilamentDensity(FilamentType type) const;
    
    // Dynamic instability (microtubules)
    void setMicrotubuleDynamics(double polymerization_rate, double depolymerization_rate,
                               double catastrophe_frequency, double rescue_frequency);
    
    // Actin dynamics
    void setActinDynamics(double polymerization_rate, double depolymerization_rate,
                         double nucleation_rate, double capping_rate);
    
    // Motor proteins
    void addMotorProtein(const std::string& motor_name, FilamentType filament,
                        double velocity, double force, double processivity);
    
    // Cargo transport
    void addCargoTransport(const std::string& cargo, const std::string& motor,
                          double loading_rate, double unloading_rate);
    
    // Mechanical properties
    void setMechanicalProperties(double elastic_modulus, double viscosity);
    double calculateMechanicalStress(double applied_force) const;
    
    void updateOrganelle(double dt, AtomSpace& atomspace, PSystem& psystem) override;
    void initializeOrganelle(AtomSpace& atomspace, PSystem& psystem) override;
    
private:
    std::map<FilamentType, double> filament_densities_;
    
    struct MicrotubuleDynamics {
        double polymerization_rate = 1.0;
        double depolymerization_rate = 0.5;
        double catastrophe_frequency = 0.1;
        double rescue_frequency = 0.05;
    } mt_dynamics_;
    
    struct ActinDynamics {
        double polymerization_rate = 2.0;
        double depolymerization_rate = 1.0;
        double nucleation_rate = 0.1;
        double capping_rate = 0.5;
    } actin_dynamics_;
    
    struct MotorProtein {
        std::string name;
        FilamentType filament;
        double velocity;     // μm/s
        double force;        // pN
        double processivity; // μm
    };
    
    std::vector<MotorProtein> motor_proteins_;
    std::map<std::string, std::pair<std::string, std::pair<double, double>>> cargo_transport_;
    
    double elastic_modulus_ = 1000.0; // Pa
    double viscosity_ = 0.1;          // Pa·s
};

/**
 * Organelle system manager
 */
class OrganelleSystem {
public:
    OrganelleSystem();
    
    // Organelle management
    void addOrganelle(std::unique_ptr<Organelle> organelle);
    Organelle* getOrganelle(const std::string& name);
    const Organelle* getOrganelle(const std::string& name) const;
    
    template<typename T>
    T* getOrganelleByType(const std::string& name) {
        auto* organelle = getOrganelle(name);
        return dynamic_cast<T*>(organelle);
    }
    
    // Initialize all organelles
    void initializeSystem(AtomSpace& atomspace, PSystem& psystem);
    
    // Simulate all organelles
    void simulateSystem(double dt, AtomSpace& atomspace, PSystem& psystem);
    
    // Inter-organelle transport
    void addInterOrganelleTransport(const std::string& source_organelle,
                                   const std::string& target_organelle,
                                   const std::string& cargo, double transport_rate);
    
    // Organelle biogenesis
    void simulateOrganelleBiogenesis(double dt, AtomSpace& atomspace);
    
    // Statistics and analysis
    std::map<std::string, double> getOrganelleVolumes() const;
    std::map<std::string, std::map<std::string, double>> getConcentrationDistribution() const;
    
private:
    std::vector<std::unique_ptr<Organelle>> organelles_;
    std::map<std::string, Organelle*> organelle_index_;
    
    struct InterOrganelleTransport {
        std::string source;
        std::string target;
        std::string cargo;
        double rate;
    };
    
    std::vector<InterOrganelleTransport> inter_organelle_transport_;
    
    void updateInterOrganelleTransport(double dt, AtomSpace& atomspace, PSystem& psystem);
};

} // namespace cellular

#endif // CELLULAR_ORGANELLES_H