# Mocks des outils externes

Faux exécutables (fastp, fastqc, multiqc, samtools, Rscript, IRFinder, SpliceLauncher.sh) qui
créent des sorties au bon nom et au bon format. Ils permettent d'exécuter les **vrais** blocs
`script:` des processus (et donc de tester l'échappement bash/Groovy, le nommage, les canaux et
les scripts de `bin/`) sans données ni outils réels :

```bash
PATH="$PWD/tests/mocks/bin:$PATH" nextflow run main.nf -profile test \
    --splicelauncher tests/mocks/SpliceLauncher/SpliceLauncher.sh
```

Ils ne valident évidemment pas les résultats biologiques.
