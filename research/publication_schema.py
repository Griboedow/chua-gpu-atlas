"""Read earlier Git publications using the current atlas field names.

Only names are translated. Decimal strings, numeric literals, timestamps and
recorded digests are retained, so forecast values can be compared exactly.
"""
import re


IDENTIFIERS = {
    'certify_chua_layer_organizer': 'certify_chua_saddle',
    'trace_chua_layer_organizer': 'estimate_chua_homoclinic_tangency',
    'chua_window_ancestry': 'chua_branch_relations',
    'refine_readme_thin_cores': 'refine_readme_thin_intervals',
    'symmetric_core_width': 'symmetric_stability_width',
    'core_beta': 'representative_beta',
    'core_sha256': 'representative_sha256',
    'core': 'representative',
    'Core': 'Representative',
    'cores': 'representatives',
    'precision_core': 'high_precision_representative',
    'advance_core': 'advance_representative',
    'new_itinerary_cores': 'new_itinerary_representatives',
    'quadratic_core_coordinate': 'normalized_representative_coordinate',
    'thin_core_samples': 'thin_interval_samples',
    'five_maximum_core_coverage': 'five_peak_stability_region_coverage',
    'first_core': 'first_stability_region',
    'second_core': 'second_stability_region',
    'max_error_over_core_width': 'max_error_over_stability_width',
    'fold_in_core_widths': 'fold_error_in_stability_widths',
    'interval_overlap_fraction_of_measured_core': 'interval_overlap_fraction',
    'predicted_midpoint_inside_measured_core': 'predicted_midpoint_inside_measured_interval',
    'raw_midpoint_in_core_units': 'raw_midpoint_in_interval_units',
    'transport_midpoint_in_core_units': 'corrected_midpoint_in_interval_units',
    'transport_fold': 'corrected_fold',
    'transport_midpoint': 'corrected_midpoint',
    'transport_midpoint_inside': 'corrected_midpoint_inside',
    'transport_fold_error_in_widths': 'corrected_fold_error_in_widths',
    'transported_midpoints_inside': 'corrected_midpoints_inside',
    'max_transport_fold_error_in_widths': 'max_corrected_fold_error_in_widths',
    'organizer': 'tangency_estimate',
    'organizer_control': 'tangency_estimate_control',
    'organizer_difference': 'tangency_estimate_difference',
    'organizer_control_difference': 'tangency_estimate_control_difference',
    'organizer_scope': 'tangency_estimate_scope',
    'ancestry': 'branch_relations',
    'local_transport': 'continue_reference_branch',
    'rejected_transport': 'rejected_continuation',
}

FILENAMES = {
    'chua-layer-prediction-organizer': 'chua-layer-prediction-tangency',
    'organizer-approximants.json': 'tangency-estimates.json',
    'transport-protocol.json': 'bias-correction-protocol.json',
    'off21-local-ancestry.json': 'off21-local-branch-relations.json',
    'offcenter-ancestry.json': 'offcenter-branch-relations.json',
    'sheet-ancestry.json': 'sheet-branch-relations.json',
    'off21-transport.json': 'off21-continuation.json',
    'thin-cores.json': 'thin-intervals.json',
}

_tokens = re.compile(r'\b(?:'+'|'.join(map(re.escape,sorted(IDENTIFIERS,key=len,reverse=True)))+r')\b')
# A recorded source-file digest describes the file at the producing revision.
# Preserve its original path as well as its digest.
_digest_path = re.compile(r'("[^"\n]+\.[a-z]+"(?=\s*:\s*"[0-9a-f]{64}"))')

LABELS = {
    'measured representative widths': 'measured stability widths',
    'representative width': 'stability width',
    'Representative width': 'Stability width',
    'shrinking representatives': 'shrinking stability intervals',
    'measured representatives at their actual width scale': 'measured stability intervals at their actual width scale',
    'numerical tangency_estimate': 'numerical tangency estimate',
    'fixed itinerary tangency_estimate': 'homoclinic tangency estimate with a fixed itinerary',
    'tangency_estimate-anchored forecast': 'tangency-based forecast',
    'improved tangency_estimate': 'improved tangency estimate',
    'coordinate_n = tangency_estimate +': 'coordinate_n = estimated_tangency_parameter +',
    'transport forecast': 'bias-corrected forecast',
    'data for transport': 'data for the bias correction',
    'transport model': 'bias-correction model',
    'cross-parameter transport': 'prediction at other parameter values',
    'locally transportable': 'approximately constant over the tested slope patch',
    'thinnest representative': 'narrowest stability interval',
    'branch branch_relations': 'branch relationships',
}


def current_names(text):
    """Translate identifiers without changing numeric data or source digests."""
    parts=_digest_path.split(text)
    for i in range(0,len(parts),2):
        for before,after in FILENAMES.items():
            parts[i]=parts[i].replace(before,after)
        parts[i]=_tokens.sub(lambda m:IDENTIFIERS[m.group()],parts[i])
        for before,after in LABELS.items():
            parts[i]=parts[i].replace(before,after)
    return ''.join(parts)
