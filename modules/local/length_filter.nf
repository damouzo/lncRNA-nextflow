// A6: Transcript length filter (>200 nt summed exon length)
// Also classifies transcripts by gffcompare class code

process TRANSCRIPT_LENGTH_FILTER {

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

    MIN_LENGTH = ${params.min_transcript_length}

    # High-confidence novel classes (kept for downstream discovery)
    KEEP_CLASSES = {'u', 'i', 'x', 'j', 'o', 'c', 'k', 'm', 'n'}
    EXCLUDE_CLASSES = {'e', 'p', 's', 'r'}  # artifact-prone, tracked separately

    def parse_attrs(attr_str):
        try:
            return dict(item.strip().split(' ', 1) for item in attr_str.rstrip(';').split('; ') if ' ' in item)
        except ValueError:
            return {}

    # Pass 1: transcript-level class code + summed exon length per transcript
    tx_class = {}
    tx_length = {}
    with open('${gffcmp_annotated_gtf}') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\\t')
            if len(parts) < 9:
                continue
            attrs = parse_attrs(parts[8])
            t_id = attrs.get('transcript_id', '').strip('"')
            if not t_id:
                continue
            if parts[2] == 'transcript':
                tx_class[t_id] = attrs.get('class_code', '').strip('"').strip()
            elif parts[2] == 'exon':
                start, end = int(parts[3]), int(parts[4])
                tx_length[t_id] = tx_length.get(t_id, 0) + (end - start + 1)

    # Transcripts retained for discovery: novel class + passing length filter
    keep = set(
        t_id for t_id in tx_class
        if tx_class[t_id] in KEEP_CLASSES and tx_length.get(t_id, 0) >= MIN_LENGTH
    )

    # Pass 2: emit full blocks (transcript/exon/CDS) for retained transcripts
    with open('${gffcmp_annotated_gtf}') as f_in, \\
         open('transcripts_length_filtered.gtf', 'w') as f_out:
        for line in f_in:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\\t')
            if len(parts) < 9:
                continue
            attrs = parse_attrs(parts[8])
            t_id = attrs.get('transcript_id', '').strip('"')
            if not t_id:
                continue
            if t_id in keep:
                f_out.write(line)

    # Class funnel — counts every class seen, not just retained ones
    total_kept = sum(1 for t_id in tx_class if t_id in keep)
    n_filtered = len(tx_class) - total_kept
    class_counts = {}
    class_kept = {}
    for t_id, code in tx_class.items():
        class_counts[code] = class_counts.get(code, 0) + 1
        if t_id in keep:
            class_kept[code] = class_kept.get(code, 0) + 1

    with open('transcripts_by_class.tsv', 'w') as out:
        out.write('class_code\\tstatus\\ttotal_transcripts\\tkept_transcripts\\n')
        for code in sorted(class_counts):
            status = ('novel' if code in KEEP_CLASSES
                      else 'artifact' if code in EXCLUDE_CLASSES
                      else 'reference')
            out.write(f"{code}\\t{status}\\t{class_counts[code]}\\t{class_kept.get(code, 0)}\\n")

    with open('filtered_stats.txt', 'w') as out:
        out.write(f"total_transcripts\\t{len(tx_class)}\\n")
        out.write(f"kept_transcripts\\t{total_kept}\\n")
        out.write(f"filtered\\t{n_filtered}\\n")
    """
}
