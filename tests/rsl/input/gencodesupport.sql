CREATE TABLE IF NOT EXISTS "gencodesupport" ("id" INTEGER NOT NULL PRIMARY KEY, "geneId" VARCHAR(255) NOT NULL, "geneName" VARCHAR(255) NOT NULL, "transcriptId" VARCHAR(255) NOT NULL, "transcriptType" VARCHAR(255) NOT NULL, "chrom" VARCHAR(255) NOT NULL, "intronStart" INTEGER NOT NULL, "intronEnd" INTEGER NOT NULL, "strand" VARCHAR(255) NOT NULL, "intronMotif" VARCHAR(255) NOT NULL, "numExprs" INTEGER NOT NULL, "numUniqueMapReads" INTEGER NOT NULL, "numMultiMapReads" INTEGER NOT NULL);
CREATE INDEX "gencodesupport_geneId" ON "gencodesupport" ("geneId");
CREATE INDEX "gencodesupport_geneName" ON "gencodesupport" ("geneName");
CREATE INDEX "gencodesupport_transcriptId" ON "gencodesupport" ("transcriptId");
CREATE INDEX "gencodesupport_transcriptType" ON "gencodesupport" ("transcriptType");
CREATE INDEX "gencodesupport_intronMotif" ON "gencodesupport" ("intronMotif");
CREATE INDEX "gencodesupport_numUniqueMapReads" ON "gencodesupport" ("numUniqueMapReads");
CREATE INDEX "gencodesupport_numMultiMapReads" ON "gencodesupport" ("numMultiMapReads");
