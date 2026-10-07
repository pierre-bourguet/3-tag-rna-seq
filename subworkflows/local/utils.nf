// Helper functions shared by the workflows

// Reference directory: --reference_dir, else <references_root>/tagseq_<te_annotation>[_<transgene>]
def resolve_reference_dir() {
    def dir = params.reference_dir ?: "${params.references_root}/${params.reference_name}"
    def manifest = file("${dir}/reference_manifest.tsv")
    if (!manifest.exists()) {
        error "No reference_manifest.tsv in ${dir}. Build it first: nextflow run main.nf --build_reference --te_annotation ${params.te_annotation}" +
              (params.transgene ? " --transgene ${params.transgene}" : '') + " -profile cbe"
    }
    return dir
}

// reference_manifest.tsv: key<TAB>value, '#' comments
def read_manifest(dir) {
    def m = [:]
    file("${dir}/reference_manifest.tsv").eachLine { line ->
        if (line.trim() && !line.startsWith('#')) {
            def f = line.split('\t', 2)
            m[f[0]] = f.size() > 1 ? f[1] : ''
        }
    }
    if (m.te_annotation && m.te_annotation != params.te_annotation && !params.reference_dir) {
        error "Reference ${dir} was built with te_annotation=${m.te_annotation}, not ${params.te_annotation}"
    }
    return m
}

// CSV with header (sample,fastq[,well,...]) or legacy headerless TSV (fastq<TAB>sample).
// Relative fastq paths are relative to the samplesheet's directory.
def samplesheet_channel(sheet) {
    def base = file(sheet).parent
    def fq = { p -> file(p.startsWith('/') ? p : "${base}/${p}", checkIfExists: true) }
    def first = file(sheet).readLines().find { it.trim() }
    def header = first.split(',')*.trim()
    if (header.containsAll(['sample', 'fastq'])) {
        return Channel.fromPath(sheet)
            .splitCsv(header: true, strip: true)
            .filter { row -> row.sample }
            .map { row -> tuple([id: row.sample, well: row.well ?: ''], fq(row.fastq)) }
    }
    return Channel.fromPath(sheet)
        .splitCsv(sep: '\t', strip: true)
        .filter { row -> row && row[0] }
        .map { row -> tuple([id: row[1], well: ''], fq(row[0])) }
}

// pipeline_info/: parameters, code version, reference manifest
def write_run_info(ref_dir) {
    def info = file("${params.outdir}/pipeline_info")
    info.mkdirs()
    def version = 'unknown'
    try {
        def p = ['git', '-C', "${projectDir}", 'describe', '--tags', '--always', '--dirty'].execute()
        p.waitFor()
        if (p.exitValue() == 0) version = p.text.trim()
    } catch (Exception e) { }
    file("${info}/run_info.tsv").text = [
        "pipeline_version\t${workflow.manifest.version}",
        "git_describe\t${version}",
        "command\t${workflow.commandLine}",
        "start\t${workflow.start}",
        "reference_dir\t${ref_dir ?: 'NA'}",
    ].join('\n') + '\n'
    file("${info}/params.json").text = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(params))
    if (ref_dir) file("${ref_dir}/reference_manifest.tsv").copyTo("${info}/reference_manifest.tsv")
}
