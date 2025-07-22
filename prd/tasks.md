# Refactoring Strategy Class

Current strategy class in `Fed-GWAS/pipeline/src/server/strategy_new.py`  is extended from FedAvg Strategy which is not good. It should be inherited from flower Strategy Class, the whole strategy class needs to be refactored. I need you to modify the strategy class to be based on flower’s Strategy class following the documentation in flower's documentation.

Links: https://flower.ai/docs/framework/how-to-implement-strategies.html

Details: 

Following the abstract class documentation of Flower Strategy Implemntation 

Implement strategies - Flower Framework 

Refactoring the current strategy class to be based on flower’s Strategy class, needs to implement the following abstract methods 

@abstractmethod 
def initialize_parameters 

@abstractmethod 
def configure_fit – specify fit instruction for each round 

@abstractmethod 
def aggregate_fit – specify server’s job at each round 

Refer to the following documentation for communication flow between client, server and strategy  

Implement strategies - Flower Framework https://flower.ai/docs/framework/how-to-implement-strategies.html
