import unittest
from scripts.main_analysis import fetch_ncbi_metadata, parse_metadata

class TestMainAnalysis(unittest.TestCase):
    def test_fetch_ncbi_metadata(self):
        query = '"H5N1"[Organism] AND "avian"[Host]'
        metadata = fetch_ncbi_metadata(query, max_results=10)
        self.assertIn("LOCUS", metadata)

    def test_parse_metadata(self):
        genbank_data = fetch_ncbi_metadata('"H5N1"[Organism] AND "avian"[Host]', 10)
        df = parse_metadata(genbank_data)
        self.assertFalse(df.empty)

if __name__ == "__main__":
    unittest.main()
