# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os

from kb_repeatmasker.repeatmasker_utils import RepeatMaskerRunner
#END_HEADER


class kb_repeatmasker:
    '''
    Module Name:
    kb_repeatmasker

    Module Description:
    A KBase module: kb_repeatmasker

Wraps RepeatMasker 4.2.4 (Smit, Hubley & Green) with the Dfam 4.0 curated
consensus library (FamDB format 3) served as KBase reference data.
Input: one or more Assemblies. Output: masked Assemblies plus a report with
repeat annotation (.out/.gff/.tbl) per input.
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.1.0"
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)
        #END_CONSTRUCTOR
        pass


    def run_repeatmasker(self, ctx, params):
        """
        :param params: instance of type "RunRepeatMaskerParams" (Parameters
           for run_repeatmasker. workspace_name - workspace for the masked
           assemblies and the report (required) assembly_refs - one or more
           assemblies to mask (required) species - NCBI Taxonomy name
           searched in Dfam, e.g. "Drosophila melanogaster", "Fungi".
           Required unless custom_library_ref is given; if both are given the
           Dfam families for the taxon and the custom library are searched
           together. custom_library_ref - Assembly of repeat consensus
           sequences, IDs ideally "name#class/family" (e.g. RepeatModeler
           output) mask_type - "soft" (lowercase, default), "hard_N" or
           "hard_X" mask_scope - "all" (default), "interspersed_only"
           (-nolow), "simple_only" (-noint) speed - "default", "quick" (-q),
           "rush" (-qq) or "slow" (-s) norna - 1 = do not mask small RNA
           (pseudo)genes (-norna) uncurated - 1 = include uncurated Dfam
           families (-uncurated) div - only mask repeats less than this %
           diverged from consensus (-div) cutoff - minimum Smith-Waterman
           score for custom-library hits (-cutoff, default 225) alignments -
           1 = also write the .align file (-a) save_masked_assembly - 1 =
           save a masked Assembly per input (default 1) output_suffix -
           appended to each input name for the masked Assembly ("_masked"))
           -> structure: parameter "workspace_name" of String, parameter
           "assembly_refs" of list of type "assembly_ref" (Reference to a
           KBaseGenomeAnnotations.Assembly (or legacy ContigSet).), parameter
           "species" of String, parameter "custom_library_ref" of type
           "assembly_ref" (Reference to a KBaseGenomeAnnotations.Assembly (or
           legacy ContigSet).), parameter "mask_type" of String, parameter
           "mask_scope" of String, parameter "speed" of String, parameter
           "norna" of type "boolean" (A boolean: 0 = false, 1 = true.),
           parameter "uncurated" of type "boolean" (A boolean: 0 = false, 1 =
           true.), parameter "div" of Double, parameter "cutoff" of Long,
           parameter "alignments" of type "boolean" (A boolean: 0 = false, 1
           = true.), parameter "save_masked_assembly" of type "boolean" (A
           boolean: 0 = false, 1 = true.), parameter "output_suffix" of String
        :returns: instance of type "RunRepeatMaskerOutput" (report_name,
           report_ref - the KBaseReport masked_assembly_refs - masked
           Assemblies, in input order repeatmasker_results - compact
           per-input summary (lengths, masked bp, classes)) -> structure:
           parameter "report_name" of String, parameter "report_ref" of
           String, parameter "masked_assembly_refs" of list of type
           "assembly_ref" (Reference to a KBaseGenomeAnnotations.Assembly (or
           legacy ContigSet).), parameter "repeatmasker_results" of list of
           unspecified object
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_repeatmasker
        output = RepeatMaskerRunner(self.config, self.callback_url, ctx).run(params)
        #END run_repeatmasker

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_repeatmasker return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
