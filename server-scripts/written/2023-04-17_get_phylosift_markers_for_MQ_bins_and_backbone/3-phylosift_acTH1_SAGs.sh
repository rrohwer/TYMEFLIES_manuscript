#!/bin/bash

/home/baker/phylosift_v1.0.1/phylosift search --isolate --besthit genomes/OID2737471675_MCM14ME106_AcTH1-A_assembled.fna && 
/home/baker/phylosift_v1.0.1/phylosift align --isolate --besthit genomes/OID2737471675_MCM14ME106_AcTH1-A_assembled.fna && 
sed 's/\.fna.*$//' PS_temp/OID2737471675_MCM14ME106_AcTH1-A_assembled.fna/alignDir/concat.updated.1.fasta > acTH1_backbone_markers/OID2737471675_MCM14ME106_AcTH1-A_assembled.faa

/home/baker/phylosift_v1.0.1/phylosift search --isolate --besthit genomes/OID2739367602__MCM14ME001_AcTH1-A_assembled.fna && 
/home/baker/phylosift_v1.0.1/phylosift align --isolate --besthit genomes/OID2739367602__MCM14ME001_AcTH1-A_assembled.fna && 
sed 's/\.fna.*$//' PS_temp/OID2739367602__MCM14ME001_AcTH1-A_assembled.fna/alignDir/concat.updated.1.fasta > acTH1_backbone_markers/OID2739367602__MCM14ME001_AcTH1-A_assembled.faa
