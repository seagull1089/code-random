default:compile

compile:
	javac Astar/src/Astar/*.java  Astar/src/test/examples/BulkTest.java -d classes
run:
	java -cp "./classes" test/examples/BulkTest MapData.txt stringPath.txt 
clean:
	rm -rf *.out classes/* 
