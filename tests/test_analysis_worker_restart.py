"""Scheduler restarts must preserve completed analysis and failed attempt records."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


WORKER = Path(__file__).resolve().parents[1]/'scripts/momentum_task.py'


class AnalysisWorkerRestart(unittest.TestCase):
    def run_case(self, exit_code=0, mutate=False):
        with tempfile.TemporaryDirectory(prefix='analysis-restart-') as directory:
            root = Path(directory)
            counter = root/'calls.txt'
            command = ['python3', '-c',
                'from pathlib import Path; import sys; p=Path(sys.argv[1]); '
                'p.write_text(p.read_text()+"call\\n" if p.exists() else "call\\n"); '
                'sys.exit(int(sys.argv[2]))', str(counter), str(exit_code)]
            task = dict(kind='analysis', command=command, output=str(root/'result'))
            tasks = root/'tasks.json'
            tasks.write_text(json.dumps([task]))

            def execute(label):
                scratch = root/label
                scratch.mkdir()
                return subprocess.run([sys.executable, str(WORKER), '--tasks', str(tasks),
                    '--index', '0'], cwd=scratch, capture_output=True, text=True)

            first = execute('first')
            self.assertEqual(first.returncode, exit_code, first.stderr)
            self.assertEqual(counter.read_text(), 'call\n')
            if mutate:
                task['command'][-1] = '9'
                tasks.write_text(json.dumps([task]))
            second = execute('second')
            self.assertEqual(second.returncode, 1 if mutate else exit_code, second.stderr)
            self.assertEqual(counter.read_text(), 'call\n')
            self.assertEqual(int((root/'result/exit_code.txt').read_text()), exit_code)
            self.assertNotEqual(json.loads((root/'result/task.json').read_text()), task) if mutate else None

    def test_success_is_reused_without_executing_again(self):
        self.run_case()

    def test_failure_is_preserved_without_executing_again(self):
        self.run_case(exit_code=7)

    def test_changed_task_cannot_reuse_existing_result(self):
        self.run_case(mutate=True)


if __name__ == '__main__':
    unittest.main()
