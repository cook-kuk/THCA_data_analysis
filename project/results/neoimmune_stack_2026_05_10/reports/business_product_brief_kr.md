# NeoImmune-Stack Business Brief KR

## Product sentence
NeoImmune-Stack is an AI operating layer that turns patient-specific mutation candidates into leakage-audited, evidence-labeled, wet-lab-ready top-N vaccine candidate queues.

## Why this is commercially sharper than another binding predictor
- Most public tools score binding/presentation; the bottleneck is immunogenicity and patient-level prioritization.
- Hospitals need a decision layer that explains why a candidate should be synthesized/tested.
- The system can absorb existing public tools without being locked to one vendor predictor.
- Local algorithms remain the differentiating science layer.

## MVP
Input: patient mutation table, HLA typing, expression, VAF/clonality, HLA LOH/APM.
Output: top-20 candidate queue, evidence rationale, rejection reasons, wet-lab test plan, claim boundary.

## First paid/partner pilot
5-10 patients with matched sequencing and HLA typing. Run all frozen public tools plus local CLEAN-Neo++ branches. Select top 20 plus controls. Validate a subset by MS/T-cell assay if available.

## Do not sell as
- validated vaccine efficacy prediction
- standalone clinical diagnostic
- LLM-generated vaccine design
