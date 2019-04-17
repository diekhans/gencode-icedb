

Commands
--------

* ``tslGetUcscRnaAligns`` - Load PSL alignments from UCSC all_mrna or all_est tables into an SQLite database. (C program)
* ``tslGetEnsemblRnaAligns`` - Fetch cDNA alignments from Ensembl and load as PSLs into an SQLite3 database.
* ``tslGbffGetProblemCases`` - Scan genbank flat files looking for known problem libraries.
* ``tslGenbankProblemCasesLoad`` - Load problem case tab file generate ``gbffGetProblemCases`` by into an SQLite3 databases.
* ``tslCollectGenbankSupport`` - Collect support for gencode annotations.  Normally run in a cluster job.
* ``tslCollectGenbankSupportMkJobs`` - Generate parasol jobs to collect TSLs for GENCODE.
* ``tslCollectGenbankSupportJob`` - Job wrapper to run ``rslGencodeCollectSupport``.
* ``tslCollectGenbankSupportFinishJobs`` - Combine ``tslCollectGenbankSupport`` job results and store in an SQLite3 table.

