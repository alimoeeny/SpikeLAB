function testPursuitIndex


%pi = PursuitIndex(MonkeyName, NeuronNumber, ClusterName, StimulusType, ti);
ti = TuningIndex('test', 001, '.c1.', 'rds', 'DT', [], 0);
pi = PursuitIndex('test', 001, '.c1.', 'rds', ti);
assertEqual(ti,-998);

ti = TuningIndex('test', 001, '.c1.', 'rds', 'DT', [], 1);
pi = PursuitIndex('test', 001, '.c1.', 'rds', ti);
assertEqual(ti,-998);

ti = TuningIndex('test', 002, '.c1.', 'rds', 'DT', [], 0);
pi = PursuitIndex('test', 002, '.c1.', 'rds', ti);
assertAlmostEqual(ti,0.8, 0.2);

ti = TuningIndex('test', 002, '.c1.', 'rds', 'DT', [], 1);
pi = PursuitIndex('test', 002, '.c1.', 'rds', ti);
assertAlmostEqual(ti,0.6, 0.2);

ti = TuningIndex('test', 003, '.c1.', 'rds', 'DT', [], 0);
pi = PursuitIndex('test', 003, '.c1.', 'rds', ti);
assertAlmostEqual(ti, -0.8, 0.2);

ti = TuningIndex('test', 003, '.c1.', 'rds', 'DT', [], 1);
pi = PursuitIndex('test', 003, '.c1.', 'rds', ti);
assertAlmostEqual(ti, -0.6, 0.2);

ti = TuningIndex('test', 004, '.c1.', 'rds', 'DT', [], 0);
pi = PursuitIndex('test', 004, '.c1.', 'rds', ti);
assertAlmostEqual(ti, 0.0, 0.2);

ti = TuningIndex('test', 004, '.c1.', 'rds', 'DT', [], 1);
pi = PursuitIndex('test', 004, '.c1.', 'rds', ti);
assertAlmostEqual(ti, 0.0, 0.1);

ti = TuningIndex('test', 005, '.c1.', 'rds', 'DT', [], 0);
pi = PursuitIndex('test', 005, '.c1.', 'rds', ti);
assertAlmostEqual(ti, 0.8, 0.1);

ti = TuningIndex('test', 005, '.c1.', 'rds', 'DT', [], 1);
pi = PursuitIndex('test', 005, '.c1.', 'rds', ti);
assertAlmostEqual(ti, 0.6, 0.2);

