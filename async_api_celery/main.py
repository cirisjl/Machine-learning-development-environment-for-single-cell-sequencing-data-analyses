import fastapi
import uvicorn
import routing

fast = fastapi.FastAPI()

fast.include_router(routing.router)


if __name__ == '__main__':
    uvicorn.run(fast, host="0.0.0.0", port=3000)

