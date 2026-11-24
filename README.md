# WaveMachine

A browser version of the **Virtual Wave Machine**, a simulation of the Shive torsional wave
apparatus, for the JHU Physics and Astronomy lecture demonstrations.

The new page (`web/`) is plain HTML and JavaScript with no libraries, build step, or server code.
Open `web/index.html` in any browser, or copy the `web/` folder to a web server.

**Instructor mode:** `index.html?instructor` adds a "Show predicted resonances" checkbox to the
energy plot. Students don't see it, since it gives away answers to the handout's Problems 3 and 4.
It's hidden, not protected.

## Credits

- **Original simulation and teaching material:** Prof. Oleg Tchernyshyov, Department of Physics and
  Astronomy, Johns Hopkins University. His Java applet (built with Easy Java Simulations) and the
  handout *Elastic waves in one dimension* (`Waves1d.pdf`) are the basis for this version. The
  physics model, parameters, intro text, and handout are his.
- **Easy Java Simulations:** Francisco Esquembre (Universidad de Murcia), with the Open Source
  Physics tools by Wolfgang Christian (Davidson College).
- **Funding for the original work:** National Science Foundation Grant No. DMR-0348679. Any
  opinions, findings, and conclusions or recommendations expressed in this material are those of
  the author and do not necessarily reflect the views of the National Science Foundation.
- **The apparatus:** invented by John N. Shive at Bell Telephone Laboratories.
- **Browser version (2026):** JHU Basic Physics Lab.

## Layout

| Folder | Contents |
|--------|----------|
| `web/` | The new simulation page (`index.html`) and its images and handout. |
| `web/instructor.html` | Instructor notes: answers to the handout problems and how the predicted resonances are calculated. Not linked from the student page. |
| `original/` | Untouched copy of the 2017 Java applet files, including the EJS source `_library/TorsionalOscillator.xml`. |
| `logs/` | Dated work log of changes and discoveries. |
| `project-rules.md` | Project goals. |

Images, PDFs, and JAR files are stored with Git LFS (see `.gitattributes`). Run `git lfs install`
once before cloning.
