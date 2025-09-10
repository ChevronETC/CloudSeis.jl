# CloudSeis.jl

| **Documentation** | **Action Statuses** |
|:---:|:---:|
| [![][docs-dev-img]][docs-dev-url] [![][docs-stable-img]][docs-stable-url] | [![][doc-build-status-img]][doc-build-status-url] [![][build-status-img]][build-status-url] [![][code-coverage-img]][code-coverage-results] |

CloudSeis.jl is a Julia library for reading and writing CloudSeis files.  The
CloudSeis data-format is designed to be similar to the JavaSeis data format while
adapting to cloud storage (e.g. Azure Blob storage).

### Debug
Make sure to update packages in Pkg if you're having dependency issues.

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://chevronetc.github.io/CloudSeis.jl/dev/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://ChevronETC.github.io/CloudSeis.jl/stable

[doc-build-status-img]: https://github.com/ChevronETC/CloudSeis.jl/workflows/Documentation/badge.svg
[doc-build-status-url]: https://github.com/ChevronETC/CloudSeis.jl/actions?query=workflow%3ADocumentation

[build-status-img]: https://github.com/ChevronETC/CloudSeis.jl/workflows/Tests/badge.svg
[build-status-url]: https://github.com/ChevronETC/CloudSeis.jl/actions?query=workflow%3A"Tests"

[code-coverage-img]: https://codecov.io/gh/ChevronETC/CloudSeis.jl/branch/master/graph/badge.svg
[code-coverage-results]: https://codecov.io/gh/ChevronETC/CloudSeis.jl