# Cellular Simulation System Documentation

## Overview

This is a comprehensive computational model for simulating human cells using OpenCog-inspired AtomSpace and membrane computing (P-systems). The framework integrates multiple scales of cellular organization from molecular interactions to organelle dynamics.

## Architecture

### Core Components

#### 1. AtomSpace (`core/atomspace.h`)
OpenCog-inspired knowledge representation system for cellular components:
- **Atoms**: Represent molecules, proteins, enzymes, organelles
- **Truth Values**: Probabilistic reasoning with strength and confidence
- **Molecular Properties**: Concentration tracking and activity levels
- **Pattern Matching**: Efficient queries for molecular interactions

#### 2. P-System (`core/psystem.h`) 
Membrane computing framework with hierarchical rules:
- **Membranes**: 13 cellular compartments (cytoplasm, nucleus, mitochondria, etc.)
- **Reaction Rules**: Mass action kinetics with stoichiometry
- **Transport Rules**: Inter-compartmental molecular transport
- **Membrane Potentials**: Goldman-Hodgkin-Katz equation implementation

### Specialized Systems

#### 3. Membrane Dynamics (`membrane/membrane_dynamics.h`)
Cell membrane structure and transport mechanisms:
- **Phospholipid Bilayer**: Asymmetric lipid distribution modeling
- **Ion Channels**: Voltage/ligand-gated channel dynamics
- **Transporters**: Michaelis-Menten kinetics for membrane transport
- **Receptors**: GPCR, RTK, and other receptor types

#### 4. Signaling Pathways (`signaling/signaling_pathways.h`)
Complete intracellular signaling networks:
- **GPCR Cascades**: G-protein coupled receptor signaling
- **Tyrosine Kinase Pathways**: RTK activation and downstream signaling
- **Second Messengers**: cAMP, IP3/DAG, calcium signaling
- **MAP Kinase Cascades**: ERK, JNK, p38 pathways

#### 5. Organelle Models (`organelles/organelles.h`)
Detailed organelle implementations:
- **Endoplasmic Reticulum**: Protein folding, calcium storage, lipid synthesis
- **Golgi Apparatus**: Protein processing and vesicle transport
- **Lysosomes**: Autophagy and hydrolytic enzyme activity
- **Peroxisomes**: ROS metabolism and fatty acid oxidation
- **Nucleus**: Transcription, RNA processing, chromatin dynamics
- **Cytoskeleton**: Filament dynamics and motor protein transport

#### 6. Mitochondria (`organelles/mitochondria.h`)
Specialized mitochondrial implementation:
- **mtDNA Operations**: Replication, transcription, maintenance
- **Electron Transport Chain**: Complexes I-IV with detailed kinetics
- **ATP Synthase**: Rotational mechanism and proton pumping
- **Membrane Dynamics**: Inner/outer membrane transport
- **Quality Control**: Mitophagy and biogenesis

#### 7. Kinetics Engine (`kinetics/kinetics_engine.h`)
Temporal dynamics and reaction kinetics:
- **Rate Laws**: Mass action, Michaelis-Menten, allosteric mechanisms
- **Integration Methods**: Euler, Runge-Kutta, adaptive stepping
- **Stochastic Simulation**: Gillespie algorithm for molecular noise
- **Multi-scale Modeling**: Fast-slow time scale separation

#### 8. Experimental Validation (`validation/experimental_validation.h`)
Validation against cell biology data:
- **Literature Database**: Standard cellular parameters and kinetics
- **Benchmark Tests**: Glucose metabolism, respiration, calcium signaling
- **Statistical Analysis**: R², RMSE, correlation analysis
- **Classic Experiments**: Warburg effect, Pasteur effect, etc.

### Integration System

#### 9. Cellular Simulation System (`cellular_simulation_system.h`)
Comprehensive integration framework:
- **Cell Types**: Generic human, hepatocyte, neuron, cardiomyocyte
- **Environmental Conditions**: Temperature, pH, oxygen, nutrients
- **Stimulus Application**: Ligands, drugs, genetic perturbations
- **Analysis Tools**: Time course analysis, flux calculations
- **Export Functions**: CSV, mathematical frameworks, reports

## Mathematical Framework

### P-System Rules
Hierarchical membrane computing rules following the mathematical framework:

```
Membrane μ: { multiset of molecules }[reaction rules]ᶜʰᵃʳᵍᵉ

Example:
Cytoplasm: { glucose(5.0), ATP(2.0), ADP(0.5) }
[glucose + ATP → glucose-6-phosphate + ADP]¹·⁰
```

### AtomSpace Representations
OpenCog-style atom representations:

```
(MOLECULE "glucose" 
  (CONCENTRATION 5.0) 
  (COMPARTMENT cytoplasm)
  (TRUTH-VALUE 1.0 0.95))

(REACTION "glycolysis_step1"
  (REACTANTS glucose ATP)
  (PRODUCTS glucose-6-phosphate ADP)
  (RATE-CONSTANT 1.5))
```

### Reaction Kinetics
Multiple rate law implementations:

1. **Mass Action**: v = k[A]ᵃ[B]ᵇ
2. **Michaelis-Menten**: v = Vₘₐₓ[S]/(Kₘ + [S])
3. **Competitive Inhibition**: v = Vₘₐₓ[S]/(Kₘ(1 + [I]/Kᵢ) + [S])
4. **Hill Equation**: v = Vₘₐₓ[S]ⁿ/(Kₐⁿ + [S]ⁿ)

### Membrane Potential
Goldman-Hodgkin-Katz equation:

```
V = (RT/F) ln((Pₖ[K⁺]ₒᵤₜ + Pₙₐ[Na⁺]ₒᵤₜ + Pᶜˡ[Cl⁻]ᵢₙ)/(Pₖ[K⁺]ᵢₙ + Pₙₐ[Na⁺]ᵢₙ + Pᶜˡ[Cl⁻]ₒᵤₜ))
```

## Usage Examples

### Basic Cell Simulation
```cpp
#include "cellular_simulation_system.h"

auto cell = CellularSimulationFactory::createGenericHumanCell();
cell->setExtracellularConditions({{"glucose", 5.0}, {"oxygen", 0.21}});
cell->runTimeCourse(60.0, 1.0);  // 60 seconds, 1s intervals
auto results = cell->getTimeCourseData();
```

### Hepatocyte Metabolism
```cpp
auto hepatocyte = CellularSimulationFactory::createHumanHepatocyte();
hepatocyte->addLigandStimulus("insulin", 1e-8, 10.0, 30.0);
hepatocyte->runSimulation();
auto fluxes = hepatocyte->getMetabolicFluxes();
```

### Experimental Validation
```cpp
ExperimentalData glucose_data;
glucose_data.experiment_name = "glucose_uptake";
glucose_data.measured_variables = {"glucose", "ATP"};
// ... set up data ...

cell->addExperimentalData(glucose_data);
double validation_score = cell->validateAgainstExperiment("glucose_uptake");
```

## Building and Testing

### Prerequisites
- C++17 compatible compiler
- Standard library with STL support

### Build Commands
```bash
cd cellular_simulation/
make all      # Build core framework
make test     # Run validation tests
make info     # Show framework components
make clean    # Clean build artifacts
```

### Test Results
The framework successfully simulates:
- ✓ Glucose metabolism (glycolysis simulation)
- ✓ ATP/ADP dynamics
- ✓ Multi-compartment transport
- ✓ Membrane potential calculations
- ✓ P-system rule execution
- ✓ AtomSpace molecular representation

## Validation Results

The system validates against experimental cell biology data:

1. **Glucose Metabolism**: Matches literature kinetics for glucose uptake and glycolysis
2. **ATP Production**: Consistent with experimental ATP/ADP ratios
3. **Membrane Transport**: Validates transporter kinetics and ion gradients
4. **Signaling Pathways**: Reproduces known pathway activation patterns
5. **Mitochondrial Function**: Matches respiratory control ratios and P/O ratios

## Deliverables Completed

✓ **Mathematical framework for membrane computing rules**
  - Hierarchical P-system implementation with 13 cellular compartments
  - Reaction and transport rules with kinetic parameters

✓ **OpenCog atomspace representations** 
  - Molecular entities with concentration and activity tracking
  - Truth values for probabilistic reasoning
  - Pattern matching for molecular interactions

✓ **Reaction rule sets for all major cellular processes**
  - Glycolysis, TCA cycle, electron transport chain
  - Signaling pathways (GPCR, RTK, MAP kinase)
  - Transport mechanisms and ion channel dynamics

✓ **Validation against experimental cell biology data**
  - Literature database with standard cellular parameters
  - Statistical validation metrics (R², RMSE, correlation)
  - Benchmark tests for major cellular processes

## Technical Implementation

### Memory Management
- Smart pointers for automatic resource management
- RAII patterns for exception safety
- Efficient STL containers for molecular storage

### Performance
- O(1) atom lookup with hash indexing
- Vectorized reaction calculations
- Adaptive time stepping for numerical stability

### Extensibility
- Template-based design for new cell types
- Plugin architecture for custom rate laws
- Modular component system

## Future Extensions

The framework provides foundation for:
- Multi-cellular tissue modeling
- Disease state simulations
- Drug discovery applications  
- Systems biology analysis
- Synthetic biology design

## References

Implementation based on established cell biology principles and computational methods from systems biology, biophysics, and membrane computing literature.