import pytest


@pytest.fixture
def tmp_path(tmp_path_factory):
    # Short names keep nested checkpoint recovery fixtures below Windows MAX_PATH.
    return tmp_path_factory.mktemp("t")
