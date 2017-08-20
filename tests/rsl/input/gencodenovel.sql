CREATE TABLE IF NOT EXISTS "gencodenovel" ("id" INTEGER NOT NULL PRIMARY KEY, "chrom" VARCHAR(255) NOT NULL, "intronStart" INTEGER NOT NULL, "intronEnd" INTEGER NOT NULL, "strand" VARCHAR(255) NOT NULL, "intronMotif" VARCHAR(255) NOT NULL, "numExprs" INTEGER NOT NULL, "numUniqueMapReads" INTEGER NOT NULL, "numMultiMapReads" INTEGER NOT NULL, "geneIds" VARCHAR(255) NOT NULL);
CREATE INDEX "gencodenovel_intronMotif" ON "gencodenovel" ("intronMotif");
CREATE INDEX "gencodenovel_numUniqueMapReads" ON "gencodenovel" ("numUniqueMapReads");
CREATE INDEX "gencodenovel_numMultiMapReads" ON "gencodenovel" ("numMultiMapReads");
CREATE INDEX "gencodenovel_geneIds" ON "gencodenovel" ("geneIds");
