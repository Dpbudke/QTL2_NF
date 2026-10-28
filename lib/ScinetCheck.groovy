// Preflight checks for runs on SCINet (Ceres) with the Slurm executor. Every entry script calls
// ScinetCheck.problems() before it starts any process and stops with the list it returns.
//
// On SCINet a project is three things at once: the Slurm account jobs are charged to, /project/<project> and
// /90daydata/<project>. A user can belong to several projects, so the run names its project once
// (--slurm_account) and every check below is made against it. Off SCINet, run with -profile docker or local,
// or pass --scinet_check false.

import java.nio.file.Path
import java.nio.file.Paths
import java.security.MessageDigest

class ScinetCheck {

    // The Nextflow release these scripts were last checked with (Ceres, October 2026).
    static final String TESTED_NEXTFLOW = '25.04'

    // usesCoordinator: main.nf and main_resume.nf, whose permutation coordinator builds paths from the pipeline
    // folder, so they must be launched there with a relative --outdir.
    static List<String> problems(Map params, def workflow, String nextflowVersion, boolean usesCoordinator = true) {
        List<String> out = []
        if (params.scinet_check == false) {
            return out
        }
        String account = params.slurm_account ?: ''
        if (!account) {
            out << 'No Slurm account. Name the SCINet project the run belongs to, for example --slurm_account do2_projects.'
            return out
        }

        // The user must belong to the project.
        String user = System.getenv('USER') ?: ''
        try {
            def proc = ['sacctmgr', '-nP', 'show', 'assoc', "user=${user}", 'format=Account'].execute()
            proc.waitFor()
            List<String> accounts = proc.text.readLines()*.trim().findAll { it }.unique()
            if (!accounts.contains(account)) {
                out << "Slurm account ${account} is not one of ${user}'s accounts (${accounts.join(', ')})."
            }
        } catch (IOException e) {
            out << 'sacctmgr is not available. Run on a SCINet login or compute node, or pass --scinet_check false off SCINet.'
        }

        // Run, work and output folders must belong to the same project.
        List<String> roots = ["/project/${account}/", "/90daydata/${account}/"]
        Path launch = workflow.launchDir.toAbsolutePath().normalize()
        Map<String, Path> folders = [
            'launch folder': launch,
            'work folder'  : workflow.workDir.toAbsolutePath().normalize(),
        ]
        String outdir = params.outdir ? params.outdir.toString() : ''
        if (outdir) {
            if (outdir.startsWith('/') && usesCoordinator) {
                out << "--outdir ${outdir} is absolute. Give it relative to the pipeline folder; the permutation coordinator builds its paths from the pipeline folder."
            } else {
                folders['output folder'] = launch.resolve(outdir).normalize()
            }
        }
        folders.each { name, path ->
            String p = path.toString() + '/'
            if (!roots.any { p.startsWith(it) }) {
                out << "The ${name} ${path} is not under ${roots.join(' or ')}."
            }
        }
        if (usesCoordinator && launch != workflow.projectDir.toAbsolutePath().normalize()) {
            out << "Launch from the pipeline folder (${workflow.projectDir}); the permutation coordinator writes under it."
        }

        // The image must be the one this version was checked with.
        String sif = params.sif ?: ''
        File image = new File(sif)
        if (!sif || !image.isFile()) {
            out << "Container image not found: ${sif}. See the README section on the container image."
        } else if (params.sif_sha256) {
            String got = sha256(image)
            if (got != params.sif_sha256) {
                out << "Container image ${sif} has SHA-256 ${got}, expected ${params.sif_sha256}."
            }
        }

        if (!nextflowVersion.startsWith(TESTED_NEXTFLOW + '.') && params.allow_untested_nextflow != true) {
            out << "Nextflow ${nextflowVersion} has not been tested with these scripts; use NXF_VER=25.04.6 (or pass --allow_untested_nextflow true)."
        }
        return out
    }

    static String sha256(File file) {
        MessageDigest md = MessageDigest.getInstance('SHA-256')
        file.withInputStream { stream ->
            byte[] buffer = new byte[8 * 1024 * 1024]
            int n
            while ((n = stream.read(buffer)) > 0) {
                md.update(buffer, 0, n)
            }
        }
        return md.digest().encodeHex().toString()
    }
}
