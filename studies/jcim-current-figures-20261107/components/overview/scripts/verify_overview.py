"""Verify freshly redrawn Figure 2 values against the approved frozen figure."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]
def numerical_view(x):
    return {'names':x['names'],
            'colours':{k:v for k,v in x['colours'].items() if k!='source'},
            'traced_example':{k:v for k,v in x['traced_example'].items() if k not in ['smiles_source','smiles_source_sha256']},
            'results':{'equal_count_grid':x['results']['equal_count_grid']['values_bacc_lo_hi'],
                       'coverage_curves':x['results']['coverage_curves']['values'],
                       'other_sources':x['results']['other_sources'],
                       'regulatory_counts':x['results']['regulatory_counts']['counts']},
            'no_refit':x['no_refit']}
actual=numerical_view(json.loads((ROOT/'outputs/overview_redraw_provenance.json').read_text()))
expected=json.loads((ROOT/'data/expected_numeric_checks.json').read_text())
assert actual==expected, 'Redrawn overview numerical/provenance values differ from the approved figure.'
assert (ROOT/'outputs/Figure_background.pdf').is_file()
assert (ROOT/'outputs/Figure_overview.pdf').is_file()
print('PASS: current Figure 2 numerical values match the approved frozen inputs; Figure 1/2 outputs exist.')
