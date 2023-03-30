# TUBRI
Data, lot selection protocol, and analyses conducted using specimens from the Tulane University Biodiversity Research Institute, in support of NSF CAREER award 2141898 entitled, "Reconstructing Parasite Abundance in River Ecosystems Over the Past Half Century"

Welcome to the GitHub repo for the TUBRI fish!

SELECTING FISH FOR DISSECTION

The entire protocol that Chelsea used to select fish for dissection can be found in lot_selection/lot_selection.R.  Here's how it worked:

I decided to focus our data collection on the lower site that Suttkus sampled on the Pearl River, where replication was best before/after 1973 and above/below the mill.  I cut out the proposed sampling on the Alabama River because there were really very few samples above the mill (see Figure 3 in the NSF CAREER proposal).  I would rather do one site really, really well instead of two shallowly.  I think this should be fine, since BACI experiments commonly involve just one site.  I also narrowed to seven fish species that were well-replicated in all four quadrants (before-control, after-control, before-impact, after-impact).

Here is how I generated list of proposed lots (shown in lot_selection/lots_suggested_for_dissection_ROUND1.csv:
- Seven species (Carpiodes velifer, Gambusia affinis, Hybognathus nuchalis, Ictalurus punctatus, Notropis atherinoides, Percina vigil, Pimephales vigilax) were selected from a list of 17 species that are common in the TUBRI collection.  Some species were eliminated because they don't occur in this stretch of the Pearl, others because they weren't well-replicated in all four quadrants (before-control, after-control, before-impact, after-impact).
- time range = 1954-2013
- I decided to focus on just one mill - the one on the lower Pearl River (i.e., the one most intensively sampled by Suttkus; see justification above).
- For each species, I aimed for n = 20 individuals per decade in each of two river stretches (n = 2, above and below the mill).
- For 7 of 9 species, I limited our request to lots with 8 or more individuals, and I ask for only 4 individual fish per lot (regardless of the number of fish in that lot).
- Ictalurus punctatus and Percina vigil are bigger fish, and there tend to be fewer per lot.  For these two species, I limited our request to lots with 4 or more individuals, and I ask for only 2 individual fish per lot.
- total number of individuals requested = 1,178 (across 332 lots)