/*
 @licstart  The following is the entire license notice for the JavaScript code in this file.

 The MIT License (MIT)

 Copyright (C) 1997-2020 by Dimitri van Heesch

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 and associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
 BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 @licend  The above is the entire license notice for the JavaScript code in this file
*/
var NAVTREE =
[
  [ "PAPRECA hybrid off-lattice kMC/MD simulator", "index.html", [
    [ "Introduction", "index.html", [
      [ "What is PAPRECA?", "index.html#intro", null ],
      [ "Authors of PAPRECA and contact details", "index.html#authors", null ],
      [ "Open-source and community guidelines", "index.html#open_source", null ],
      [ "Acknowledgments and citations", "index.html#citations", null ],
      [ "Bibliography", "index.html#bibliography", null ]
    ] ],
    [ "Download and Installation", "md_installation.html", [
      [ "PAPRECA repository and download", "md_installation.html#download", [
        [ "PAPRECA autobuild", "md_installation.html#autobuild", null ]
      ] ],
      [ "Read before you install MANUALLY: Important notes and prerequisites", "md_installation.html#beforeInstall", null ],
      [ "Building PAPRECA manually", "md_installation.html#build", [
        [ "CMake", "md_installation.html#cmake", null ],
        [ "Traditional Make", "md_installation.html#traditionalMake", null ]
      ] ]
    ] ],
    [ "Essential information (See before you start)", "md_essentials.html", [
      [ "PAPRECA flowchart", "md_essentials.html#flowchart", null ],
      [ "Supported Predefined events", "md_essentials.html#supEvents", null ],
      [ "Domain Decomposition", "md_essentials.html#decomposition", null ],
      [ "Random Numbers and Repeatability", "md_essentials.html#rnums", null ],
      [ "PAPRECA and LAMMPS input Files", "md_essentials.html#inputFiles", null ],
      [ "Unsupported Features and illegal commands in LAMMPS input files.", "md_essentials.html#unsupfeaut", null ],
      [ "Dealing with large log.lammps files", "md_essentials.html#loglammps", null ],
      [ "Running a PAPRECA simulation", "md_essentials.html#running", null ],
      [ "fix papreca", "md_essentials.html#fixpapreca", null ],
      [ "PAPRECA units", "md_essentials.html#units", null ],
      [ "Bibliography", "md_essentials.html#essentials_bibliography", null ]
    ] ],
    [ "PAPRECA Commands", "md_commands.html", [
      [ "fix papreca command", "md_commands.html#FIX_papreca", [
        [ "Syntax", "md_commands.html#FIX_papreca_syntax", null ],
        [ "Description", "md_commands.html#FIX_papreca_description", null ]
      ] ],
      [ "kMC_steps command", "md_commands.html#KMC_steps", [
        [ "Syntax", "md_commands.html#KMCsteps_syntax", null ],
        [ "Example(s)", "md_commands.html#KMCsteps_examples", null ],
        [ "Description", "md_commands.html#KMCsteps_description", null ]
      ] ],
      [ "KMC_per_MD command", "md_commands.html#KMC_per_MD", [
        [ "Syntax", "md_commands.html#KMC_per_MD_syntax", null ],
        [ "Example(s)", "md_commands.html#KMC_per_MD_examples", null ],
        [ "Description", "md_commands.html#KMC_per_MD_description", null ]
      ] ],
      [ "random_seed command", "md_commands.html#random_seed", [
        [ "Syntax", "md_commands.html#random_seed_syntax", null ],
        [ "Example(s)", "md_commands.html#random_seed_examples", null ],
        [ "Description", "md_commands.html#random_seed_description", null ],
        [ "Bibliography", "md_commands.html#random_seed_bibliography", null ]
      ] ],
      [ "sigmas_options command", "md_commands.html#sigoptions", [
        [ "Syntax", "md_commands.html#sigoptions_syntax", null ],
        [ "Example(s)", "md_commands.html#sigoptions_examples", null ],
        [ "Description", "md_commands.html#sigoptions_description", null ],
        [ "Bibliography", "md_commands.html#sigoptions_bibliography", null ]
      ] ],
      [ "init_sigma command", "md_commands.html#insigma", [
        [ "Syntax", "md_commands.html#insigma_syntax", null ],
        [ "Example(s)", "md_commands.html#insigma_examples", null ],
        [ "Description", "md_commands.html#insigma_description", null ],
        [ "Bibliography", "md_commands.html#insigma_bibliography", null ]
      ] ],
      [ "fluid_atomtypes command", "md_commands.html#flutypes", [
        [ "Syntax", "md_commands.html#flutypes_syntax", null ],
        [ "Example(s)", "md_commands.html#flutypes_examples", null ],
        [ "Description", "md_commands.html#flutypes_description", null ]
      ] ],
      [ "frozen_atomtypes command", "md_commands.html#frotypes", [
        [ "Syntax", "md_commands.html#frotypes_syntax", null ],
        [ "Example(s)", "md_commands.html#frotypes_examples", null ],
        [ "Description", "md_commands.html#frotypes_description", null ],
        [ "Default", "md_commands.html#frotypes_defaults", null ]
      ] ],
      [ "time_end command", "md_commands.html#time_end", [
        [ "Syntax", "md_commands.html#time_end_syntax", null ],
        [ "Example(s)", "md_commands.html#time_end_examples", null ],
        [ "Description", "md_commands.html#time_end_description", null ],
        [ "Default", "md_commands.html#time_end_default", null ]
      ] ],
      [ "height_calculation command", "md_commands.html#height_calculation", [
        [ "Syntax", "md_commands.html#height_calculation_syntax", null ],
        [ "Example(s)", "md_commands.html#height_calculation_examples", null ],
        [ "Description", "md_commands.html#height_calculation_description", null ],
        [ "Default", "md_commands.html#height_calculation_default", null ]
      ] ],
      [ "desorption command", "md_commands.html#desorb", [
        [ "Syntax", "md_commands.html#desorb_syntax", null ],
        [ "Example(s)", "md_commands.html#desorb_examples", null ],
        [ "Description", "md_commands.html#desorb_description", null ],
        [ "Default", "md_commands.html#desorb_default", null ]
      ] ],
      [ "minimize_prior command", "md_commands.html#minprior", [
        [ "Syntax", "md_commands.html#minprior_syntax", null ],
        [ "Example(s)", "md_commands.html#minprior_examples", null ],
        [ "Description", "md_commands.html#minprior_description", null ],
        [ "Default", "md_commands.html#minprior_default", null ],
        [ "Bibliography", "md_commands.html#minprior_bibliography", null ]
      ] ],
      [ "minimize_after command", "md_commands.html#minafter", [
        [ "Syntax", "md_commands.html#minafter_syntax", null ],
        [ "Example(s)", "md_commands.html#minafter_examples", null ],
        [ "Description", "md_commands.html#minafter_description", null ],
        [ "Default", "md_commands.html#minafter_default", null ],
        [ "Bibliography", "md_commands.html#minafter_bibliography", null ]
      ] ],
      [ "trajectory_duration command", "md_commands.html#trajdur", [
        [ "Syntax", "md_commands.html#trajdur_syntax", null ],
        [ "Example(s)", "md_commands.html#trajdur_examples", null ],
        [ "Description", "md_commands.html#trajdur_description", null ],
        [ "Default", "md_commands.html#trajdur_default", null ]
      ] ],
      [ "create_BondBreak command", "md_commands.html#createBreak", [
        [ "Syntax", "md_commands.html#createBreak_syntax", null ],
        [ "Example(s)", "md_commands.html#createBreak_examples", null ],
        [ "Description", "md_commands.html#createBreak_description", null ]
      ] ],
      [ "create_BondForm command", "md_commands.html#createForm", [
        [ "Syntax", "md_commands.html#createForm_syntax", null ],
        [ "Example(s)", "md_commands.html#createForm_examples", null ],
        [ "Description", "md_commands.html#createForm_description", null ]
      ] ],
      [ "species_maxbonds command", "md_commands.html#maxbonds", [
        [ "Syntax", "md_commands.html#maxbonds_syntax", null ],
        [ "Example(s)", "md_commands.html#maxbonds_examples", null ],
        [ "Description", "md_commands.html#maxbonds_description", null ],
        [ "Default", "md_commands.html#maxbonds_default", null ]
      ] ],
      [ "species_maxbondtypes command", "md_commands.html#maxbondtypes", [
        [ "Syntax", "md_commands.html#maxbondtypes_syntax", null ],
        [ "Example(s)", "md_commands.html#maxbondtypes_examples", null ],
        [ "Description", "md_commands.html#maxbondtypes_description", null ],
        [ "Default", "md_commands.html#maxbondtypes_default", null ]
      ] ],
      [ "create_Deposition command", "md_commands.html#createDepo", [
        [ "Syntax", "md_commands.html#createDepo_syntax", null ],
        [ "Example(s)", "md_commands.html#createDepo_examples", null ],
        [ "Description", "md_commands.html#createDepo_description", null ],
        [ "Default", "md_commands.html#createDepo_defaults", null ],
        [ "Bibliography", "md_commands.html#createDepo_bibliography", null ]
      ] ],
      [ "random_depovecs command", "md_commands.html#depovecs", [
        [ "Syntax", "md_commands.html#depovecs_syntax", null ],
        [ "Example(s)", "md_commands.html#depovecs_examples", null ],
        [ "Description", "md_commands.html#depovevecs_description", null ],
        [ "Default", "md_commands.html#depovecs_default", null ]
      ] ],
      [ "depoheights command", "md_commands.html#depoheights", [
        [ "Syntax", "md_commands.html#depoheights_syntax", null ],
        [ "Example(s)", "md_commands.html#depoheights_examples", null ],
        [ "Description", "md_commands.html#depoheights_description", null ],
        [ "Default", "md_commands.html#depoheights_default", null ]
      ] ],
      [ "create_DiffusionHop command", "md_commands.html#createDiff", [
        [ "Syntax", "md_commands.html#createDiff_syntax", null ],
        [ "Example(s)", "md_commands.html#createDiff_examples", null ],
        [ "Description", "md_commands.html#createDiff_description", null ],
        [ "Bibliography", "md_commands.html#createDiff_bibliography", null ]
      ] ],
      [ "random_diffvecs command", "md_commands.html#diffvecs", [
        [ "Syntax", "md_commands.html#diffvecs_syntax", null ],
        [ "Example(s)", "md_commands.html#diffvecs_examples", null ],
        [ "Description", "md_commands.html#diffvecs_description", null ],
        [ "Default", "md_commands.html#diffvecs_default", null ]
      ] ],
      [ "create_MonoatomicDesorption command", "md_commands.html#createMonodes", [
        [ "Syntax", "md_commands.html#createMonodes_syntax", null ],
        [ "Example(s)", "md_commands.html#createMonodes_examples", null ],
        [ "Description", "md_commands.html#createMonodes_description", null ]
      ] ],
      [ "export_HeightVtime command", "md_commands.html#heightvtime", [
        [ "Syntax", "md_commands.html#heightvtime_syntax", null ],
        [ "Example(s)", "md_commands.html#heightVtime_examples", null ],
        [ "Description", "md_commands.html#heightVtime_description", null ],
        [ "Default", "md_commands.html#heightVtime_default", null ]
      ] ],
      [ "export_SurfaceCoverage command", "md_commands.html#coverage", [
        [ "Syntax", "md_commands.html#coverage_syntax", null ],
        [ "Example(s)", "md_commands.html#coverage_examples", null ],
        [ "Description", "md_commands.html#coverage_description", null ],
        [ "Default", "md_commands.html#coverage_default", null ]
      ] ],
      [ "export_ElementalDistributions command", "md_commands.html#Edistributions", [
        [ "Syntax", "md_commands.html#Edistributions_syntax", null ],
        [ "Example(s)", "md_commands.html#Edistributions_examples", null ],
        [ "Description", "md_commands.html#Edistributions_description", null ],
        [ "Default", "md_commands.html#Edistributions_default", null ]
      ] ],
      [ "export_ExecutionTimes command", "md_commands.html#execution", [
        [ "Syntax", "md_commands.html#execution_syntax", null ],
        [ "Example(s)", "md_commands.html#execution_examples", null ],
        [ "Description", "md_commands.html#execution_description", null ],
        [ "Default", "md_commands.html#execution_default", null ]
      ] ],
      [ "restart_freq command", "md_commands.html#restart", [
        [ "Syntax", "md_commands.html#restart_syntax", null ],
        [ "Example(s)", "md_commands.html#restart_examples", null ],
        [ "Description", "md_commands.html#restart_description", null ],
        [ "Default", "md_commands.html#restart_default", null ]
      ] ]
    ] ],
    [ "Example Applications", "md_examples.html", [
      [ "Your first PAPRECA simulation: Brownian Diffusion (random walk)", "md_examples.html#brownian", [
        [ "LAMMPS and PAPRECA input files", "md_examples.html#brownian_INPUT", null ],
        [ "Execution", "md_examples.html#brownian_run", null ],
        [ "Results", "md_examples.html#brownian_results", null ],
        [ "Bibliography", "md_examples.html#brownian_bibliography", null ]
      ] ],
      [ "Monoatomic adsorption and desorption", "md_examples.html#adsorption", [
        [ "LAMMPS and PAPRECA input files", "md_examples.html#adsorption_INPUT", null ],
        [ "Execution", "md_examples.html#adsorption_run", null ],
        [ "Results", "md_examples.html#adsorption_results", null ],
        [ "Bibliography", "md_examples.html#adsorption_bibliography", null ]
      ] ],
      [ "Film growth from the thermal decomposition of tricresyl phosphate (TCP) molecules on a Fe110 surface.", "md_examples.html#phosphates", [
        [ "LAMMPS and PAPRECA input files", "md_examples.html#phosphates_INPUT", null ],
        [ "Execution", "md_examples.html#phosphates_run", null ],
        [ "Results", "md_examples.html#phosphates_results", null ],
        [ "Bibliography", "md_examples.html#phosphates_bibliography", null ]
      ] ],
      [ "Organic solvents", "md_examples.html#solvents", [
        [ "Execution", "md_examples.html#solvents_run", null ]
      ] ]
    ] ],
    [ "Namespaces", "namespaces.html", [
      [ "Namespace List", "namespaces.html", "namespaces_dup" ],
      [ "Namespace Members", "namespacemembers.html", [
        [ "All", "namespacemembers.html", null ],
        [ "Functions", "namespacemembers_func.html", null ],
        [ "Typedefs", "namespacemembers_type.html", null ]
      ] ]
    ] ],
    [ "Classes", "annotated.html", [
      [ "Class List", "annotated.html", "annotated_dup" ],
      [ "Class Index", "classes.html", null ],
      [ "Class Hierarchy", "hierarchy.html", "hierarchy" ],
      [ "Class Members", "functions.html", [
        [ "All", "functions.html", "functions_dup" ],
        [ "Functions", "functions_func.html", "functions_func" ],
        [ "Variables", "functions_vars.html", null ],
        [ "Related Functions", "functions_rela.html", null ]
      ] ]
    ] ],
    [ "Files", "files.html", [
      [ "File List", "files.html", "files_dup" ],
      [ "File Members", "globals.html", [
        [ "All", "globals.html", null ],
        [ "Functions", "globals_func.html", null ],
        [ "Macros", "globals_defs.html", null ]
      ] ]
    ] ]
  ] ]
];

var NAVTREEINDEX =
[
"annotated.html",
"classPAPRECA_1_1PaprecaConfig.html#a9bb2ab226ce712074a229abcd718e5e2",
"event__detect_8cpp.html#a7ddc46c695cdcb6604bb7ccc15cfe6fe",
"input__file_8h.html#af5284300eb1d6e5b966950617cac9b29",
"namespacePAPRECA.html#a187a68bba2930ef8b13a802dca683cc2",
"utilities_8h.html#a90c8c7ad8a281214e4abb1c9fc5d4f31"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';