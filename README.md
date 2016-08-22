# The Statistical Scribbles blog

A graphical and code-based walk through statistical concepts that I am currently teaching.


## How was it set-up?
The blog is served as a page from my main website at [openpencil.github.io](http://openpencil.github.io/).
To make Github recognize the contents of the `statistical_scribbles` repository as a page on the main website,
the repository should have only the `gh-pages` branch and no `master` branch.
- Create the gh-pages branch
`git checkout -b gh-pages`
`git push origin gh-pages`
- In the Github GUI, switch the "Default Branch" under repository settings to `gh-pages`.
- Switch to `gh-pages` branch, and delete the local master branch
`git branch -d master`
- Delete the remote master branch with:
`git push origin :master`

If there are no errors, the contents of the repository (a jekyll blog) will be served as a page from the main website.

## TODO
- Document specific Jekyll settings
- Figure out how to eliminate fluff from design
- Make the initial posts

## Sources
- The initial source for this whittled-down repository came from: [Jekyll with tiny-feet](http://rpiai.com/jekyll-tinyfeet/)
