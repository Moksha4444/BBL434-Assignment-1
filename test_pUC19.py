# test_pUC19.py
"""
Test suite for plasmid constructor
Validates restriction site removal and sequence integrity
"""

from plasmid_builder import PlasmidAssembler


def test_pUC19_restriction_site_removal():
    """Verify that restriction sites are properly removed from backbone"""
    
    # Construct plasmid using the assembler
    plasmid_sequence = PlasmidAssembler.construct_plasmid(
        "pUC19.fa", 
        "Design_pUC19.txt"
    )
    
    # Test 1: EcoRI site must be removed from backbone
    assert "GAATTC" not in plasmid_sequence, "❌ FAIL: EcoRI site still present in backbone!"
    
    # Test 2: Validate DNA sequence integrity
    gc_count = plasmid_sequence.count("G") + plasmid_sequence.count("C")
    assert gc_count > 0, "❌ FAIL: Invalid DNA sequence - no G or C bases"
    
    # Test 3: Ensure sequence is not empty
    assert len(plasmid_sequence) > 0, "❌ FAIL: Empty plasmid sequence"
    
    # Test 4: Check for valid DNA bases only
    valid_bases = set("ATGC")
    sequence_bases = set(plasmid_sequence)
    assert sequence_bases.issubset(valid_bases), f"❌ FAIL: Invalid bases found: {sequence_bases - valid_bases}"
    
    print("✓ All tests passed!")
    print(f"✓ EcoRI successfully removed from backbone")
    print(f"✓ Plasmid length: {len(plasmid_sequence)} bp")
    print(f"✓ GC content: {(gc_count / len(plasmid_sequence) * 100):.1f}%")


if __name__ == "__main__":
    test_pUC19_restriction_site_removal()