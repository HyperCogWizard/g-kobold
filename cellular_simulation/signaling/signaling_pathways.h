#ifndef CELLULAR_SIGNALING_PATHWAYS_H
#define CELLULAR_SIGNALING_PATHWAYS_H

#include "../core/atomspace.h"
#include "../core/psystem.h"
#include <functional>
#include <map>
#include <vector>

namespace cellular {

/**
 * G-protein types
 */
enum class GProteinType {
    GS,     // Stimulatory (activates adenylyl cyclase)
    GQ_G11, // Activates phospholipase C
    GI_GO,  // Inhibitory (inhibits adenylyl cyclase)
    G12_G13 // Activates Rho signaling
};

/**
 * Second messenger types
 */
enum class SecondMessengerType {
    CAMP,           // Cyclic adenosine monophosphate
    CGMP,           // Cyclic guanosine monophosphate
    IP3,            // Inositol 1,4,5-trisphosphate
    DAG,            // Diacylglycerol
    CALCIUM,        // Calcium ions
    NITRIC_OXIDE    // Nitric oxide
};

/**
 * Protein kinase types
 */
enum class KinaseType {
    PKA,        // Protein kinase A (cAMP-dependent)
    PKC,        // Protein kinase C (DAG/Ca2+ dependent)
    PKG,        // Protein kinase G (cGMP-dependent)
    CAMKII,     // Ca2+/calmodulin-dependent protein kinase II
    ERK,        // Extracellular signal-regulated kinase
    JNK,        // c-Jun N-terminal kinase
    P38,        // p38 MAP kinase
    AKT,        // AKT/PKB kinase
    TYROSINE_KINASE // Receptor tyrosine kinase
};

/**
 * G-protein coupled receptor cascade
 */
class GPCRCascade {
public:
    GPCRCascade(const std::string& receptor_name, GProteinType g_protein);
    
    const std::string& getReceptorName() const { return receptor_name_; }
    GProteinType getGProteinType() const { return g_protein_type_; }
    
    // Receptor binding and activation
    void setLigandBinding(double kd, double max_binding);
    double calculateReceptorActivation(double ligand_concentration) const;
    
    // G-protein activation
    void setGProteinProperties(double activation_rate, double gtp_hydrolysis_rate);
    double calculateGProteinActivity(double receptor_activity, double gtp_concentration) const;
    
    // Downstream effector activation
    void addEffector(const std::string& effector_name, double activation_constant);
    std::map<std::string, double> calculateEffectorActivation(double g_protein_activity) const;
    
    // Update signaling state
    void updateSignaling(double dt, AtomSpace& atomspace, PSystem& psystem);
    
private:
    std::string receptor_name_;
    GProteinType g_protein_type_;
    double kd_ = 1e-9;              // Ligand binding affinity
    double max_binding_ = 1000.0;   // Maximum binding sites
    double g_activation_rate_ = 1.0; // G-protein activation rate
    double gtp_hydrolysis_rate_ = 0.1; // GTP hydrolysis rate
    std::map<std::string, double> effectors_; // Downstream effectors
};

/**
 * Tyrosine kinase signaling pathway
 */
class TyrosineKinasePathway {
public:
    TyrosineKinasePathway(const std::string& receptor_name);
    
    const std::string& getReceptorName() const { return receptor_name_; }
    
    // Receptor dimerization and autophosphorylation
    void setDimerizationProperties(double kd_dimer, double phosphorylation_rate);
    double calculateReceptorDimerization(double ligand_concentration) const;
    double calculateAutophosphorylation(double dimer_concentration) const;
    
    // SH2/SH3 domain interactions
    void addAdapterProtein(const std::string& protein, double binding_affinity);
    std::map<std::string, double> calculateAdapterBinding(double phospho_receptor) const;
    
    // Downstream kinase cascades
    void addDownstreamKinase(KinaseType kinase, double activation_threshold);
    std::map<KinaseType, double> calculateKinaseActivation(
        const std::map<std::string, double>& adapter_proteins) const;
    
    void updateSignaling(double dt, AtomSpace& atomspace, PSystem& psystem);
    
private:
    std::string receptor_name_;
    double kd_dimer_ = 1e-6;
    double phosphorylation_rate_ = 1.0;
    std::map<std::string, double> adapter_proteins_;
    std::map<KinaseType, double> downstream_kinases_;
};

/**
 * Second messenger system
 */
class SecondMessengerSystem {
public:
    SecondMessengerSystem();
    
    // cAMP signaling
    void setCyclicAMPProperties(double adenylyl_cyclase_rate, double phosphodiesterase_rate);
    double calculatecAMPLevel(double gs_protein_activity, double current_camp) const;
    
    // IP3/DAG signaling
    void setIP3DAGProperties(double plc_rate, double ip3_degradation_rate);
    std::pair<double, double> calculateIP3DAG(double gq_protein_activity, double pip2_level) const;
    
    // Calcium signaling
    void setCalciumProperties(double release_rate, double uptake_rate, double extrusion_rate);
    double calculateCalciumRelease(double ip3_level, double er_calcium) const;
    double calculateCalciumDynamics(double current_ca, double release_rate, double dt) const;
    
    // Protein kinase activation
    double calculatePKAActivity(double camp_level) const;
    double calculatePKCActivity(double dag_level, double calcium_level) const;
    
    void updateSecondMessengers(double dt, AtomSpace& atomspace, PSystem& psystem);
    
private:
    double adenylyl_cyclase_rate_ = 1.0;
    double phosphodiesterase_rate_ = 0.5;
    double plc_rate_ = 1.0;
    double ip3_degradation_rate_ = 2.0;
    double ca_release_rate_ = 5.0;
    double ca_uptake_rate_ = 2.0;
    double ca_extrusion_rate_ = 1.0;
};

/**
 * Calcium signaling network
 */
class CalciumSignalingNetwork {
public:
    CalciumSignalingNetwork();
    
    // Calcium sources
    void setExtracellularCalcium(double concentration) { extracellular_ca_ = concentration; }
    void setERCalciumStores(double concentration) { er_ca_stores_ = concentration; }
    
    // Calcium channels
    void addCalciumChannel(const std::string& channel_name, double conductance, double threshold);
    
    // Calcium pumps and exchangers
    void addCalciumPump(const std::string& pump_name, double max_rate, double km);
    
    // Calcium buffering
    void addCalciumBuffer(const std::string& buffer_name, double kd, double capacity);
    
    // Calcium-dependent processes
    void addCalciumTarget(const std::string& target_name, double kd, double hill_coeff);
    
    // Simulation
    double calculateCalciumInflux(double membrane_potential, double driving_force) const;
    double calculateCalciumEfflux(double cytosolic_calcium) const;
    double calculateBuffering(double free_calcium) const;
    std::map<std::string, double> calculateCalciumTargets(double calcium_level) const;
    
    void updateCalciumDynamics(double dt, AtomSpace& atomspace, PSystem& psystem);
    
private:
    double extracellular_ca_ = 2.0; // mM
    double er_ca_stores_ = 0.5;     // mM
    std::map<std::string, std::pair<double, double>> calcium_channels_; // conductance, threshold
    std::map<std::string, std::pair<double, double>> calcium_pumps_;    // max_rate, km
    std::map<std::string, std::pair<double, double>> calcium_buffers_;  // kd, capacity
    std::map<std::string, std::pair<double, double>> calcium_targets_;  // kd, hill_coeff
};

/**
 * MAP kinase cascades
 */
class MAPKinaseCascade {
public:
    MAPKinaseCascade(const std::string& pathway_name);
    
    const std::string& getPathwayName() const { return pathway_name_; }
    
    // Three-tier kinase cascade: MAPKKK -> MAPKK -> MAPK
    void setKinaseTier(int tier, const std::string& kinase_name, 
                      double activation_rate, double deactivation_rate);
    
    // Scaffold proteins
    void addScaffoldProtein(const std::string& scaffold_name, double binding_affinity);
    
    // Upstream activators
    void addUpstreamActivator(const std::string& activator_name, double activation_strength);
    
    // Downstream targets
    void addDownstreamTarget(const std::string& target_name, double phosphorylation_rate);
    
    // Cascade simulation
    std::vector<double> calculateCascadeActivation(
        const std::map<std::string, double>& upstream_signals, double dt) const;
    
    std::map<std::string, double> calculateTargetPhosphorylation(
        const std::vector<double>& kinase_activities) const;
    
    void updateCascade(double dt, AtomSpace& atomspace, PSystem& psystem);
    
private:
    std::string pathway_name_;
    struct KinaseTier {
        std::string name;
        double activation_rate = 1.0;
        double deactivation_rate = 0.5;
        double activity = 0.0;
    };
    std::vector<KinaseTier> kinase_tiers_; // MAPKKK, MAPKK, MAPK
    std::map<std::string, double> scaffold_proteins_;
    std::map<std::string, double> upstream_activators_;
    std::map<std::string, double> downstream_targets_;
};

/**
 * Comprehensive intracellular signaling system
 */
class IntracellularSignaling {
public:
    IntracellularSignaling();
    
    // Pathway management
    void addGPCRCascade(const GPCRCascade& cascade);
    void addTyrosineKinasePathway(const TyrosineKinasePathway& pathway);
    void addMAPKinaseCascade(const MAPKinaseCascade& cascade);
    
    // Signaling components
    SecondMessengerSystem& getSecondMessengerSystem() { return second_messengers_; }
    CalciumSignalingNetwork& getCalciumNetwork() { return calcium_network_; }
    
    // Cross-talk and integration
    void addCrossTalk(const std::string& source_pathway, const std::string& target_pathway,
                     double interaction_strength);
    
    // Temporal dynamics
    void addFeedbackLoop(const std::string& pathway, const std::string& target,
                        double feedback_strength, bool positive_feedback);
    
    // Simulation
    void simulateSignaling(double dt, AtomSpace& atomspace, PSystem& psystem);
    
    // Analysis
    std::map<std::string, double> getPathwayActivities() const;
    std::map<SecondMessengerType, double> getSecondMessengerLevels() const;
    
private:
    std::vector<GPCRCascade> gpcr_cascades_;
    std::vector<TyrosineKinasePathway> rtk_pathways_;
    std::vector<MAPKinaseCascade> mapk_cascades_;
    SecondMessengerSystem second_messengers_;
    CalciumSignalingNetwork calcium_network_;
    
    // Cross-talk matrix
    std::map<std::pair<std::string, std::string>, double> cross_talk_;
    
    // Feedback loops
    std::map<std::pair<std::string, std::string>, std::pair<double, bool>> feedback_loops_;
    
    void updateCrossTalk(AtomSpace& atomspace, PSystem& psystem);
    void updateFeedbackLoops(AtomSpace& atomspace, PSystem& psystem);
};

} // namespace cellular

#endif // CELLULAR_SIGNALING_PATHWAYS_H