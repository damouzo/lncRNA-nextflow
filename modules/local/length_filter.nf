// A6: Transcript length filter (>200 nt summed exon length)
// Also classifies transcripts by gffcompare class code

process TRANSCRIPT_LENGTH_FILTER {
    label 'small_task'

    input:
    path gffcmp_annotated_gtf
    path class_summary

    output:
    path "transcripts_length_filtered.gtf", emit: filtered_gtf
    path "transcripts_by_class.tsv", emit: class_table
    path "filtered_stats.txt", emit: stats

    script:
    """
    #!/usr/bin/env python3
    import sys

    MIN_LENGTH = ${params.min_transcript_length}

    # High-confidence novel classes
    KEEP_CLASSES = {'u', 'i', 'x', 'j', 'o', 'c', 'k', 'm', 'n'}
    EXCLUDE_CLASSES = {'e', 'p', 's', 'r'}  # artifact-prone, tracked separately

    transcripts = {}
    with open('${gffcmp_annotated_gtf}') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\\t')
            if len(parts) < 9:
                continue
            if parts[2] != 'exon':
                continue
            attrs = dict(item.strip().split(' ', 1) for item in parts[8].rstrip(';').split('; ') if ' ' in item)
            t_id = attrs.get('transcript_id', '').strip('"')
            if not t_id:
                continue
            start, end = int(parts[3]), int(parts[4])
            exon_len = end - start + 1
            transcripts[t_id] = transcripts.get(t_id, 0) + exon_len

    n_filtered = 0
    n_kept = 0
    with open('${gffcmp_annotated_gtf}') as f_in, \\
         open('transcripts_length_filtered.gtf', 'w') as f_out:
        for line in f_in:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\\t')
            if len(parts) < 9:
                continue
            if parts[2] == 'transcript':
                attrs = dict(item.strip().split(' ', 1) for item in parts[8].rstrip(';').split('; ') if ' ' in item)
                t_id = attrs.get('transcript_id', '').strip('"')
                t_len = transcripts.get(t_id, 0)
                if t_len >= MIN_LENGTH:
                    n_kept += 1
                    f_out.write(line)
                else:
                    n_filtered += 1

    with open('filtered_stats.txt', 'w') as out:
        out.write(f"total_transcripts\\t{n_kept + n_filtered}\\n")
        out.write(f"kept_transcripts\\t{n_kept}\\n")
        out.write(f"filtered_length\\t{n_filtered}\\n")
    """
}
