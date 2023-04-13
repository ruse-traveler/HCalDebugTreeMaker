#!/bin/bash
# -----------------------------------------------------------------------------
# 'RunHCalTreeMaker.sh'
# Derek Anderson
# 04.13.2023
#
# A simple script to run the
# JBarrelHCalTreeMaker JANA
# plugin
# -----------------------------------------------------------------------------

# i/o files
input="../forTestingClusterTree.e10th7080n1evt100pim.d13m4y2023.edm4hep.root"
podio="forTestingClusterTree.e10th7080n1evt100pim.d13m4y2023.podio.root"
output="forTestingClusterTree.e10th7080n1evt100pim.d13m4y2023.plugin.root"

# output collections from EICrecon
collections="HcalBarrelRecHits,HcalBarrelClusters,HcalBarrelIslandProtoClusters,HcalBarrelTruthClusters,HcalBarrelTruthProtoClusters,GeneratedParticles"

eicrecon -Pplugins=JBarrelHCalTreeMaker -Ppodio:output_include_collections=$collections -Peicrecon:LogLevel=debug -Ppodio:output_file=$podio -Phistsfile=$output $input

# end -------------------------------------------------------------------------
