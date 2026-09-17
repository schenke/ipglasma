########################################################
#
#    Doxygen developer-reference documentation.
#
########################################################

find_package(Doxygen OPTIONAL_COMPONENTS dot)

if (DOXYGEN_FOUND)
    # In the following, CMake variables are set that look unused here but are
    # substituted into Doxyfile.in by configure_file() below.
    set(DOXYGEN_EXTRACT_ALL "YES")

    if (Doxygen_dot_FOUND)
        set(HAVE_DOT "YES")
        get_target_property(DOT_EXECUTABLE Doxygen::dot IMPORTED_LOCATION)
        get_filename_component(DOT_PATH "${DOT_EXECUTABLE}" DIRECTORY)
    else ()
        set(HAVE_DOT "NO")
        set(DOT_PATH "")
        message(ATTENTION "Graphviz 'dot' not found -- diagrams won't be generated in the documentation.")
    endif ()

    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
                    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

    add_custom_target(doc
                      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
                      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                      COMMENT "Building Doxygen developer documentation in ${CMAKE_CURRENT_BINARY_DIR}/html"
                      VERBATIM)

    # A second configuration, with EXTRACT_ALL off, makes Doxygen itself
    # enumerate every documentation gap (WARN_IF_UNDOCUMENTED only fires for
    # entities that would otherwise not appear in the docs at all). Used by
    # the three targets below to check/enforce completeness; none of them
    # write out real HTML, so also skip dot for speed.
    set(DOXYGEN_EXTRACT_ALL "NO")
    set(HAVE_DOT "NO")
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
                    ${CMAKE_CURRENT_BINARY_DIR}/DoxyfileOnlyDocumented @ONLY)

    # Matches only the "missing documentation" family of Doxygen warnings
    # (member/parameter/return-value not documented), not incidental parse
    # warnings unrelated to documentation coverage (e.g. Markdown files using
    # GitHub-style $...$ math, which Doxygen doesn't understand and reports
    # as "unknown command" -- see README.md's note on why it stays that way).
    set(UNDOCUMENTED_PATTERN "(is|are) not documented")

    find_program(GREP grep)
    if (GREP)
        add_custom_target(undocumented
                          ${DOXYGEN_EXECUTABLE}
                          ${CMAKE_CURRENT_BINARY_DIR}/DoxyfileOnlyDocumented
                          2>&1
                          |
                          ${GREP}
                          -iE
                          "${UNDOCUMENTED_PATTERN}"
                          ||
                          ${CMAKE_COMMAND} -E echo "No undocumented entities found."
                          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                          COMMENT "Listing every Doxygen 'not documented' warning"
                          VERBATIM)

        add_custom_target(undocumented_count
                          ${DOXYGEN_EXECUTABLE}
                          ${CMAKE_CURRENT_BINARY_DIR}/DoxyfileOnlyDocumented
                          2>&1
                          |
                          ${GREP}
                          -icE
                          "${UNDOCUMENTED_PATTERN}"
                          |
                          xargs
                          -I{}
                          ${CMAKE_COMMAND} -E echo "There are {} Doxygen warnings about missing documentation."
                          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                          COMMENT "Counting Doxygen 'not documented' warnings"
                          VERBATIM)

        # The CI-enforced completeness gate: `!` negates grep's exit status,
        # so this whole pipeline fails (non-zero exit) as soon as grep finds
        # at least one matching line, i.e. at least one undocumented entity.
        add_custom_target(undocumented_test
                          !
                          ${DOXYGEN_EXECUTABLE}
                          ${CMAKE_CURRENT_BINARY_DIR}/DoxyfileOnlyDocumented
                          2>&1
                          |
                          ${GREP}
                          -qiE
                          "${UNDOCUMENTED_PATTERN}"
                          WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                          COMMENT "Testing documentation completeness for CI: fails if there is any warning"
                          VERBATIM)
    endif ()

    mark_as_advanced(GREP)
else ()
    message(ATTENTION "Doxygen not found -- the 'doc' target won't be available.")
endif ()
