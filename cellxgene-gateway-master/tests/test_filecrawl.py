import unittest
from collections import defaultdict
from unittest.mock import MagicMock, patch

from cellxgene_gateway.filecrawl import (
    render_item,
    render_item_source,
    render_item_tree,
)
from cellxgene_gateway.items.file.fileitem import FileItem
from cellxgene_gateway.items.file.fileitem_source import FileItemSource
from cellxgene_gateway.items.item import ItemTree, ItemType

source = FileItemSource("/tmp")


class TestRenderEntry(unittest.TestCase):
    def test_GIVEN_path_both_slash_THEN_view_has_single_slash(self):
        entry = FileItem(subpath="/somepath/", name="entry", type=ItemType.h5ad)
        rendered = render_item(entry, source)
        self.assertIn("view/somepath/entry/'", rendered)

    def test_GIVEN_path_starts_slash_THEN_view_has_single_slash(self):
        entry = FileItem(subpath="/somepath", name="entry", type=ItemType.h5ad)
        rendered = render_item(entry, source)
        self.assertIn("view/somepath/entry/'", rendered)

    def test_GIVEN_path_ends_slash_THEN_view_has_single_slash(self):
        entry = FileItem(subpath="somepath/", name="entry", type=ItemType.h5ad)
        rendered = render_item(entry, source)
        self.assertIn("view/somepath/entry/'", rendered)

    def test_GIVEN_path_no_slash_THEN_view_has_single_slash(self):
        entry = FileItem(subpath="somepath", name="entry", type=ItemType.h5ad)
        rendered = render_item(entry, source)
        self.assertIn("view/somepath/entry/'", rendered)


class TestRenderItemSource(unittest.TestCase):
    @patch("cellxgene_gateway.items.file.fileitem_source.FileItemSource")
    def test_GIVEN_some_filter_THEN_includes_filterpart_in_heading(self, item_source):
        item_source.name = "FakeSource"
        item_source.list_items.return_value = ItemTree("rootdir", [], [])
        rendered = render_item_source(item_source, "some_filter")
        self.assertEqual(
            rendered,
            "<h6><a href='/filecrawl.html?source=FakeSource'>FakeSource</a>:some_filter</h6><li><a href='/filecrawl/rootdir?source=FakeSource'>rootdir</a><ul></ul></li>",
        )


class TestRenderItemTree(unittest.TestCase):
    def setUp(self):
        from cellxgene_gateway.gateway import app

        self.app = app
        self.app_context = self.app.test_request_context()
        self.app_context.push()

    @patch("cellxgene_gateway.items.file.fileitem_source.FileItemSource")
    def test_GIVEN_deep_nested_dirs_THEN_includes_dirs_in_output(self, item_source):
        item_source.name = "FakeSource"
        item_source.get_annotations_subpath = lambda _: "FakeAnnotations"
        file_item = FileItem(
            subpath="foo/bar/baz", name="file.h5ad", type=ItemType.h5ad
        )
        item_tree = ItemTree("foo/bar/baz", [file_item], [])
        rendered = render_item_tree(item_tree, item_source)
        self.assertEqual(
            rendered,
            "<li><a href='/filecrawl/foo/bar/baz?source=FakeSource'>baz</a><ul>"
            "<li> <a href='/source/FakeSource/view/foo/bar/baz/file.h5ad/'>file.h5ad</a>"
            "  | annotations: <a class='new' href='/source/FakeSource/view/FakeAnnotations'>new</a>"
            "</li></ul></li>",
        )

    @patch(
        "os.listdir",
        side_effect=lambda parent: defaultdict(
            list, {"tmp": ["foo"], "tmp/foo": ["bar"]}
        )[parent],
    )
    @patch("os.path.exists", return_value=True)
    def test_GIVEN_dirs_without_h5ad_THEN_excludes_dirs_in_output(
        self, listdir, exists
    ):
        # Directories:
        # - tmp
        #   - foo
        #     - bar (no h5ad files)
        item_source = FileItemSource("tmp", name="local")
        item_tree = item_source.list_items("foo")
        rendered = render_item_tree(item_tree, item_source)
        self.assertEqual(
            rendered,
            "<li><a href='/filecrawl/foo?source=local'>foo</a><ul></ul></li>",
        )
