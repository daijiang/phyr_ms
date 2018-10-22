# Manuscript for the `phyr` package

[![Travis build status](https://travis-ci.org/daijiang/phyr_ms.svg?branch=master)](https://travis-ci.org/daijiang/phyr_ms)

After pushing to Github, Travis will automatically try to build the new `ms.html` file for us (it may take up to 10 mins after the push). If you changes do not affect `ms.Rmd`, append `[skip travis]` in your git commit messge to avoid rebuilding. The `ms.html` file can be viewed at here: <https://htmlpreview.github.io/?https://github.com/daijiang/phyr_ms/blob/master/ms.html>.

You can also update `ms.html` in your computer. To render the Rmd file into HTML file, go to Terminal, `cd` change directory to the folder of this manuscript, then type `make`. Then you browser should open the file for you after it is done. You can render to MS Word with `make ms.docx` too if you want (I see no reason for this at this stage).

Citations in `ms.Rmd`:

1. Go to [Google Scholar](https://scholar.google.com); click the top left, and select `Settings`; then in the `Bibliography manager` section, check `show links to import citations into BibTeX`, then Save it. This step only need to be done once.
2. Back to the Google Scholar homepage and search the paper you want to cite. E.g. "Generalized linear mixed models for phylogenetic analyses of community structure", under the item, click `Import into BibTeX`; then copy the whole new page (`Cmd + A` then `Cmd + C`).

    ```
    @article{ives2011generalized,
      title={Generalized linear mixed models for phylogenetic analyses of community structure},
      author={Ives, Anthony R and Helmus, Matthew R},
      journal={Ecological Monographs},
      volume={81},
      number={3},
      pages={511--525},
      year={2011},
      publisher={Wiley Online Library}
    }
    ```

3. Open the `ref.bib` in this folder, and paste the citation information there.
4. To cite the paper in the manuscript, copy its key (e.g. `ives2011generalized` here), and inset `[@ives2011generalized]` to where you want to cite it. For multiple papers, use `[@key1; @key2]`.
5. That's it.

