function testNeurClus

[MonkeyName, NeuronNumber, ClusterName] = NeurClus({'dae062'});
assertTrue(strcmp(MonkeyName, 'dae'));
assertEqual(NeuronNumber, 62);
assertTrue(strcmp(ClusterName, '.c1.'));


[MonkeyName, NeuronNumber, ClusterName] = NeurClus({'dae0622'});
assertTrue(strcmp(MonkeyName, 'dae'));
assertEqual(NeuronNumber, 62);
assertTrue(strcmp(ClusterName, '.c2.'));

[MonkeyName, NeuronNumber, ClusterName] = NeurClus({'ic6222'});
assertTrue(strcmp(MonkeyName, 'icarus'));
assertEqual(NeuronNumber, 622);
assertTrue(strcmp(ClusterName, '.c2.'));

[MonkeyName, NeuronNumber, ClusterName] = NeurClus({'ic623'});
assertTrue(strcmp(MonkeyName, 'icarus'));
assertEqual(NeuronNumber, 623);
assertTrue(strcmp(ClusterName, '.c1.'));

[MonkeyName, NeuronNumber, ClusterName] = NeurClus({'ic458'});
assertTrue(strcmp(MonkeyName, 'icarus'));
assertEqual(NeuronNumber, 458);
assertTrue(strcmp(ClusterName, '.c1.'));

[MonkeyName, NeuronNumber, ClusterName] = NeurClus({'tst001'});
assertTrue(strcmp(MonkeyName, 'test'));
assertEqual(NeuronNumber, 1);
assertTrue(strcmp(ClusterName, '.c1.'));


[MonkeyName, NeuronNumber, ClusterName] = NeurClus({'tst123'});
assertTrue(strcmp(MonkeyName, 'test'));
assertEqual(NeuronNumber, 123);
assertTrue(strcmp(ClusterName, '.c1.'));

[MonkeyName, NeuronNumber, ClusterName] = NeurClus({'tst1232'});
assertTrue(strcmp(MonkeyName, 'test'));
assertEqual(NeuronNumber, 123);
assertTrue(strcmp(ClusterName, '.c2.'));
