import awkward as ak

from hzupsilonphoton.builders import (
#    build_boson,
#    build_bosons_combination,
#    build_dimuons,
    build_good_muons,
    build_good_photons,
#    build_mu_1,
#    build_mu_2,
#    build_photon,
#    build_upsilon,
    build_probe_muon,
    build_tag_muon,
    build_probe_photon,
    build_TrigObj,
)
from hzupsilonphoton.feed_forward import (
    FeedForwardSequence,
    FilterSequence,
    ObjectSequence,
    WeightSequence,
)
from hzupsilonphoton.filters import (
    lumisection_filter,
#    mass_selection_filter,
#    signal_selection_filter,
    trigger_filter,
    n_muons_filter,
    n_probe_muon_filter,
    n_probe_photon_filter,
    n_tag_muon_filter,
)
from hzupsilonphoton.weighters import (
    generator_weight,
    l1prefr_weights,
#    muon_id_weight,
#    muon_iso_weight,
#    photon_electron_veto_weight,
#    photon_id_weight,
    pileup_weight,
)

forward_events = FeedForwardSequence("base_sequence")
forward_events.register_sequence(FilterSequence("lumisection", lumisection_filter))
forward_events.register_sequence(FilterSequence("trigger", trigger_filter))

forward_events.register_sequence(WeightSequence("pileup", pileup_weight))
forward_events.register_sequence(WeightSequence("generator", generator_weight))
forward_events.register_sequence(WeightSequence("l1_prefiring", l1prefr_weights))

forward_events.register_sequence(ObjectSequence("good_muons", build_good_muons))
forward_events.register_sequence(ObjectSequence("good_photons", build_good_photons))

forward_events.register_sequence(
    FilterSequence("n_muons", lambda evts: ak.num(evts.events.good_muons) >= 2)
)
forward_events.register_sequence(
    FilterSequence("n_photons", lambda evts: ak.num(evts.events.good_photons) >= 1)
)
#forward_events.register_sequence(ObjectSequence("dimuons", build_dimuons))
#forward_events.register_sequence(
#    FilterSequence("n_dimuons", lambda evts: ak.num(evts.events.dimuons) >= 1)
#)


forward_events.register_sequence(ObjectSequence("probe_photon", build_probe_photon))
forward_events.register_sequence(
    FilterSequence("n_probe_photon", lambda evts: ak.num(evts.events.probe_photon) >= 1)
)

#forward_events.register_sequence(WeightSequence("pileup", pileup_weight))
#forward_events.register_sequence(WeightSequence("generator", generator_weight))
#forward_events.register_sequence(WeightSequence("l1_prefiring", l1prefr_weights))


forward_events.register_sequence(ObjectSequence("probe_muon", build_probe_muon))
forward_events.register_sequence(ObjectSequence("tag_muon", build_tag_muon))
forward_events.register_sequence(ObjectSequence("TrigObj", build_TrigObj))

forward_events.register_sequence(
    FilterSequence("n_probe_muon", lambda evts: ak.num(evts.events.probe_muon) >= 2)
)

forward_events.register_sequence(
    FilterSequence("n_tag_muon", lambda evts: ak.num(evts.events.tag_muon) >= 1)
)

