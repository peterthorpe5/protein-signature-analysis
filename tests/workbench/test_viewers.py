"""Scientific display contracts for cropped alignments and labelled phylograms."""

import pytest

from protein_signature_app.workbench_viewers import alignment_window, tree_figure
from protein_signatures.errors import InputValidationError


def test_alignment_window_is_an_explicit_crop():
    """Display bounds preserve original gaps and one-based alignment coordinates."""
    result = alignment_window(text=">first\nAC-D\n>second\nACED\n", start=2, width=2)
    assert "C-" in result and "CE" in result
    assert result.count("\n") == 1


@pytest.mark.parametrize(
    "text,start",
    [
        ("", 1),
        ("ACD", 1),
        (">a\nAC\n>b\nA", 1),
        (">a\nAC", 3),
        (">a\nAC", 0),
        (">a\nAC\n>a\nAC", 1),
    ],
)
def test_alignment_view_rejects_invalid_matrices_and_windows(text, start):
    """Malformed matrices and out-of-range coordinates cannot produce misleading previews."""
    with pytest.raises(InputValidationError):
        alignment_window(text=text, start=start)


def test_phylogram_retains_exact_tip_mapping_and_branch_scale():
    """Plot labels must link safe tool identifiers back to actual source proteins."""
    figure = tree_figure(
        newick="(a:0.1,(b:0.2,c:0.3):0.1);", labels={"a": "P1|long", "b": "NR_2.1", "c": "plant:C"}
    )
    assert set(figure.data[1].text) == {"P1|long", "NR_2.1", "plant:C"}
    assert "substitutions" in figure.layout.xaxis.title.text
    zero = tree_figure(newick="(a:0,b:0,c:0);", labels={"a": "A", "b": "B", "c": "C"})
    assert "all inferred branch lengths were zero" in zero.layout.xaxis.title.text


@pytest.mark.parametrize(
    "newick,labels",
    [
        ("(a,b,c)", {"a": "A"}),
        ("(a:0.1,b:0.2,c:0.3);", {"a": "A"}),
        ("(a:-1,b:0.1,c:0.2);", {"a": "A", "b": "B", "c": "C"}),
    ],
)
def test_phylogram_rejects_incomplete_mappings_and_invalid_lengths(newick, labels):
    """A display cannot silently omit tips or treat invalid branches as valid evidence."""
    with pytest.raises(InputValidationError):
        tree_figure(newick=newick, labels=labels)
