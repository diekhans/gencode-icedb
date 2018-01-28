SELECT cmp.chrom, cmp.txStart, cmp.txEnd, att.geneId, att.transcriptId, att.geneName, att.geneType, tsl.level
FROM wgEncodeGencodeAttrsV27 AS att,
     wgEncodeGencodeTranscriptionSupportLevelV27 AS tsl,
     wgEncodeGencodeCompV27 AS cmp
WHERE att.transcriptId = tsl.transcriptId
      AND cmp.name = att.transcriptId
ORDER BY chrom, txStart, strand, geneId;
