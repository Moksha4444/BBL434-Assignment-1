#!/usr/bin/env python3
#Universal plasmid constructor

import sys
import math
import os
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class RestrictionEnzymeDatabase:
    #sequences from the markers.tab file given as the input form the professor
    
    ENZYMES = {
        "EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT",
        "PstI": "CTGCAG", "KpnI": "GGTACC", "SacI": "GAGCTC",
        "SalI": "GTCGAC", "SphI": "GCATGC", "XbaI": "TCTAGA",
        "NotI": "GCGGCCGC", "SmaI": "CCCGGG", "BsaI": "GGTCTC",
        "BbsI": "GAAGAC", "BsmBI": "CGTCTC"
    }
    
    #recognition site sequence for the given enzyme
    @classmethod
    def get_site(cls, enzyme_name):
        return cls.ENZYMES.get(enzyme_name, None)
    
    #all enzymes names
    @classmethod
    def all_enzymes(cls):
        return cls.ENZYMES.keys()


class SequenceAnalyzer:
    
    @staticmethod
    def calculate_nucleotide_frequencies(sequence):
        #Calculate background nucleotide probabilities
        length = len(sequence)
        frequencies = {
            "A": sequence.count("A") / length,
            "C": sequence.count("C") / length,
            "G": sequence.count("G") / length,
            "T": sequence.count("T") / length
        }
        return frequencies
    
    #filter out low complexity kmers
    @staticmethod
    def assess_sequence_complexity(kmer):
        unique_bases = len(set(kmer))
        return unique_bases > 2
    
    @staticmethod
    def calculate_cumulative_gc_skew(dna_sequence):
        skew_values = [0]
        running_total = 0
        
        for nucleotide in dna_sequence:
            if nucleotide == "G":
                running_total += 1
            elif nucleotide == "C":
                running_total -= 1
            skew_values.append(running_total)
        
        return skew_values
    
    #ori is approximated at the minimum of cumulative GC skew
    @staticmethod
    def locate_replication_origin_position(dna_sequence):
        skew_profile = SequenceAnalyzer.calculate_cumulative_gc_skew(dna_sequence)
        minimum_position = skew_profile.index(min(skew_profile))
        return minimum_position

#position weight matrix as taught in class
class MotifPWM:
    
    def __init__(self, motif_instances, motif_length):
        self.length = motif_length
        self.matrix = self._construct_matrix(motif_instances)
    
    def _construct_matrix(self, instances):
        #Build PWM from motif instances with pseudocounts
        matrix = []
        
        for position in range(self.length):
            nucleotide_counts = {"A": 1, "C": 1, "G": 1, "T": 1}
            
            for motif in instances:
                nucleotide_counts[motif[position]] += 1
            
            total_counts = sum(nucleotide_counts.values())
            position_probs = {base: count / total_counts 
                            for base, count in nucleotide_counts.items()}
            matrix.append(position_probs)
        
        return matrix
    
    def calculate_sequence_score(self, sequence, background_probs):
        total_score = 0
        
        for idx in range(len(sequence) - self.length + 1):
            subsequence = sequence[idx:idx + self.length]
            log_likelihood = 0
            
            for pos, base in enumerate(subsequence):
                motif_prob = self.matrix[pos][base]
                bg_prob = background_probs[base]
                log_likelihood += math.log(motif_prob / bg_prob)
            
            total_score += log_likelihood
        
        return total_score

#origin finder class
class OriginFinder:
    
    def __init__(self, k_range=(6, 15)):
        self.min_k = k_range[0]
        self.max_k = k_range[1]
    
    #background neucleotide model for likelihood calculation
    def discover_motifs(self, sequence, kmer_size):
        bg_model = SequenceAnalyzer.calculate_nucleotide_frequencies(sequence)
        optimal_score = float("-inf")
        optimal_motifs = None
        
        for start_pos in range(len(sequence) - kmer_size + 1):
            seed_motif = sequence[start_pos:start_pos + kmer_size]
            
            if not SequenceAnalyzer.assess_sequence_complexity(seed_motif):
                continue
            
            candidate_motifs = [seed_motif]
            
            for scan_pos in range(len(sequence) - kmer_size + 1):
                candidate = sequence[scan_pos:scan_pos + kmer_size]
                if SequenceAnalyzer.assess_sequence_complexity(candidate):
                    candidate_motifs.append(candidate)
            
            pwm_model = MotifPWM(candidate_motifs, kmer_size)
            current_score = pwm_model.calculate_sequence_score(sequence, bg_model)
            
            if current_score > optimal_score:
                optimal_score = current_score
                optimal_motifs = candidate_motifs
        
        return optimal_motifs, optimal_score
 
#For finding the optimal k mer length and finding the D as well
    def select_optimal_kmer_length(self, sequence):
        evaluation_results = []
        
        for k_value in range(self.min_k, self.max_k + 1):
            motif_set, score = self.discover_motifs(sequence, k_value)
            evaluation_results.append((k_value, score, motif_set))
        
        best_result = max(evaluation_results, key=lambda x: x[1])
        return best_result
    
    def identify_origin_region(self, fasta_path, window_size=500):
        #Complete ORI identification pipeline
        genome_record = SeqIO.read(fasta_path, "fasta")
        genome_sequence = str(genome_record.seq).upper()
        
        #Locate ORI center using GC skew
        origin_center = SequenceAnalyzer.locate_replication_origin_position(genome_sequence)
        
        # Extract region around ORI
        region_start = max(0, origin_center - window_size // 2)
        region_end = min(len(genome_sequence), origin_center + window_size // 2)
        origin_region = genome_sequence[region_start:region_end]
        
        # Analyze motifs in region
        best_k, likelihood_score, discovered_motifs = self.select_optimal_kmer_length(origin_region)
        
        return {
            "center_position": origin_center,
            "region_bounds": (region_start, region_end),
            "extracted_sequence": origin_region,
            "optimal_kmer": best_k,
            "likelihood": likelihood_score,
            "motif_collection": discovered_motifs[:10]
        }


class PlasmidAssembler:
    #Constructs custom plasmid from components
    
    @staticmethod
    def load_genetic_marker(marker_identifier):
        #Retrieve marker sequence from file
        marker_file = f"markers/{marker_identifier}.fa"
        
        if not os.path.exists(marker_file):
            print(f"⚠ Missing marker file: {marker_identifier}")
            return ""
        
        marker_record = SeqIO.read(marker_file, "fasta")
        return str(marker_record.seq).upper()
    
    @staticmethod
    def parse_plasmid_design(design_filepath):
        #Extract components from design specification
        cloning_sites = []
        resistance_genes = []
        reporter_genes = []
        
        with open(design_filepath, 'r') as design_file:
            for line in design_file:
                line_content = line.strip()
                if not line_content:
                    continue
                
                component_type, component_name = line_content.split(",")
                component_type = component_type.strip().lower()
                component_name = component_name.strip()
                
                if "site" in component_type:
                    cloning_sites.append(component_name)
                elif "gene" in component_type:
                    resistance_genes.append(component_name)
                else:
                    reporter_genes.append(component_name)
        
        return cloning_sites, resistance_genes, reporter_genes
    
    @staticmethod
    def sanitize_backbone_sequence(sequence, enzyme_list):
        #Remove unwanted restriction sites from functional regions
        cleaned_sequence = sequence
        
        for enzyme in enzyme_list:
            recognition_site = RestrictionEnzymeDatabase.get_site(enzyme)
            if recognition_site:
                cleaned_sequence = cleaned_sequence.replace(recognition_site, "")
        
        return cleaned_sequence
    
    @classmethod
    def construct_plasmid(cls, input_genome, design_spec):
        #Main assembly pipeline
        
        #Phase 1: Identify and extract origin
        ori_finder = OriginFinder()
        origin_data = ori_finder.identify_origin_region(input_genome)
        
        origin_sequence = origin_data["extracted_sequence"]
        start_coord, end_coord = origin_data["region_bounds"]
        
        #Display origin discovery results
        print("\n═══ ORIGIN IDENTIFICATION ═══")
        print(f"→ Central position: {origin_data['center_position']}")
        print(f"→ Extracted region: {start_coord} to {end_coord}")
        print(f"→ Optimal k-mer: {origin_data['optimal_kmer']}")
        
        #Phase 2: Parse design requirements
        mcs_enzymes, antibiotic_markers, screening_markers = cls.parse_plasmid_design(design_spec)
        
        #Phase 3: Assemble backbone
        plasmid_backbone = origin_sequence
        
        #Add resistance markers
        for resistance_gene in antibiotic_markers:
            plasmid_backbone += cls.load_genetic_marker(resistance_gene)
        
        #Add screening markers
        for screening_gene in screening_markers:
            plasmid_backbone += cls.load_genetic_marker(screening_gene)
        
        #Phase 4: Clean backbone (remove conflicting sites)
        plasmid_backbone = cls.sanitize_backbone_sequence(
            plasmid_backbone, 
            RestrictionEnzymeDatabase.all_enzymes()
        )
        
        #Phase 5: Add MCS (must be preserved)
        for enzyme in mcs_enzymes:
            site_sequence = RestrictionEnzymeDatabase.get_site(enzyme)
            if not site_sequence:
                print(f"⚠ Unknown enzyme: {enzyme}")
                continue
            plasmid_backbone += site_sequence
        
        return plasmid_backbone


def main():
    #Application entry point
    
    if len(sys.argv) != 3:
        print("Usage: python plasmid_constructor.py <genome.fa> <design.txt>")
        sys.exit(1)
    
    genome_file = sys.argv[1]
    design_file = sys.argv[2]
    
    #Construct plasmid
    final_plasmid = PlasmidAssembler.construct_plasmid(genome_file, design_file)
    
    #Generate output
    output_record = SeqRecord(
        Seq(final_plasmid),
        id="Engineered_Plasmid",
        description="Custom plasmid with identified ORI and selected markers"
    )
    
    output_filename = "Output.fa"
    SeqIO.write(output_record, output_filename, "fasta")
    
    print(f"✓ Plasmid sequence saved to {output_filename}")
    print(f"✓ Total length: {len(final_plasmid)} bp\n")


if __name__ == "__main__":
    main()