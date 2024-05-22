* use declare variant
  (https://www.openmp.org/wp-content/uploads/openmp-examples-5.2.1.pdf page
  494) for info/warn etc errors which do nothing
* Declare variant candidates:
  - Messages::Warn (and other printers)
  - DefUtils::GetCurrentElement (threading focused branching and threadprivate variable)
  - DefUtils::GetElementNodes / DefUtils::GetElementNodesVec
* Pessimistic about GPU-ability:
  - DefUtils::GetElementNOFDOFS (mostly branching, very little arithmetic computations). If this is kept on cpu then there is no point GPUing GetCurrentElement

# 2024-04-26

* Adding some `omp target enter data directives will speed execution lot
  - Without any `data` mapping directives it is SLOW
  - With the `enter data` we get 4-0.5 execution time of cpu

# 2024-05-08

* TODO: 
  - [ ] Test cpu perf on 035, 035x2 and 015 meshes and compare to gpu perf. Make note here: 
    - [ ] Pretty ok when looking at profiler (~300us per color for 015 mesh)
  - [ ] Glueing to global crs matrix
  - [ ] p-elements
    - ElementInfoVec:iin lippu "IsReferenceElement". Tällöin ei lasketa mitään ylimääräistä (ltogmap yms), vaan ainoastaan suunnistukset CPU:lla.
    - Lasketaan kantafunktiot kerran ja asetetaan paikoilleen. Sitten GPU:lla mäpätään ltogmap yms jutut.
  - [o] Monta mpi rankkia?
    - Dedikoidaan yksi rank / node, siirretään assemblydata sinne ja assembloidaan ja palautetaan rankeille:
      MPI_COMM_SPLIT_TYPE(old_comm, split_type, key, info, new_comm). split_type = MPI_COMM_TYPE_SHARED 
    - [o] Voiko koota assemblydatan rankeittain ja kopioida vain pointterin yhdelle rankille joka ajaa omp target regionin?
    - Ei järkeä:
      - Verkko pitää värittää jotta sen voi koota säikeistetysti
      - paitsi jos laittaa pari rankkia per gpu (esim 4) jolloin väritys ei veisi liian pitkään
  - [o] mpi rank per gpu ? 
  - Rocalution: CSR col-pointer pitää olla järjestetty (nouseva)
  - `void rocalution::LocalMatrix::MultiColoring(int &num_colors, int **size_colors, LocalVector<int> *permutation) const`
