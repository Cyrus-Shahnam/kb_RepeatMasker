/*
A KBase module: kb_repeatmasker

Wraps RepeatMasker 4.2.4 (Smit, Hubley & Green) with the Dfam 4.0 curated
consensus library (FamDB format 3) served as KBase reference data.
Input: one or more Assemblies. Output: masked Assemblies plus a report with
repeat annotation (.out/.gff/.tbl) per input.
*/
module kb_repeatmasker {

    /* A boolean: 0 = false, 1 = true. */
    typedef int boolean;

    /* Reference to a KBaseGenomeAnnotations.Assembly (or legacy ContigSet). */
    typedef string assembly_ref;

    /*
    Parameters for run_repeatmasker.

    workspace_name - workspace for the masked assemblies and the report (required)
    assembly_refs - one or more assemblies to mask (required)
    species - NCBI Taxonomy name searched in Dfam, e.g. "Drosophila melanogaster",
        "Fungi". Required unless custom_library_ref is given; if both are given the
        Dfam families for the taxon and the custom library are searched together.
    custom_library_ref - Assembly of repeat consensus sequences, IDs ideally
        "name#class/family" (e.g. RepeatModeler output)
    mask_type - "soft" (lowercase, default), "hard_N" or "hard_X"
    mask_scope - "all" (default), "interspersed_only" (-nolow), "simple_only" (-noint)
    speed - "default", "quick" (-q), "rush" (-qq) or "slow" (-s)
    norna - 1 = do not mask small RNA (pseudo)genes (-norna)
    uncurated - 1 = include uncurated Dfam families (-uncurated)
    div - only mask repeats less than this % diverged from consensus (-div)
    cutoff - minimum Smith-Waterman score for custom-library hits (-cutoff, default 225)
    alignments - 1 = also write the .align file (-a)
    save_masked_assembly - 1 = save a masked Assembly per input (default 1)
    output_suffix - appended to each input name for the masked Assembly ("_masked")
    */
    typedef structure {
        string workspace_name;
        list<assembly_ref> assembly_refs;
        string species;
        assembly_ref custom_library_ref;
        string mask_type;
        string mask_scope;
        string speed;
        boolean norna;
        boolean uncurated;
        float div;
        int cutoff;
        boolean alignments;
        boolean save_masked_assembly;
        string output_suffix;
    } RunRepeatMaskerParams;

    /*
    report_name, report_ref - the KBaseReport
    masked_assembly_refs - masked Assemblies, in input order
    repeatmasker_results - compact per-input summary (lengths, masked bp, classes)
    */
    typedef structure {
        string report_name;
        string report_ref;
        list<assembly_ref> masked_assembly_refs;
        list<UnspecifiedObject> repeatmasker_results;
    } RunRepeatMaskerOutput;

    funcdef run_repeatmasker(RunRepeatMaskerParams params)
        returns (RunRepeatMaskerOutput output) authentication required;
};
